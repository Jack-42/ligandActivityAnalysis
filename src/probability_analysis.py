"""
@author Jack Ringer
Date: 7/14/2025
Description:
This code is related to gather_similarity_values.py, and attempts to answer the question:
'If two ligands have similarity > 0.85, what is the probability they share at least one cluster?'.
Credit to Claude for helping write some of the code here.
"""

import argparse
import os
import sys
from itertools import combinations

import numpy as np
import polars as pl
from tqdm import tqdm

from utils.constants import get_ligand2cluster_fpath
from utils.io import load_from_pkl
from utils.ltm import check_size, get_id2idx_map, get_ltm_idx


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute probability values of interest for multiple clusters",
        epilog="",
    )
    parser.add_argument(
        "--cluster_dir",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input directory where clustered ligands and targets were saved to by assign_family_clusters.py",
    )
    parser.add_argument(
        "--similarity_npy_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output .npy file containing lower-triangular similarity matrix",
    )
    parser.add_argument(
        "--similarity_id_pkl_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output .pkl file containing compound IDs (in order used by similarity_npy_file)",
    )
    parser.add_argument(
        "--prob_analysis_dir",
        type=str,
        required=True,
        help="Output directory where .pkl and .tsv files containing computed probabilities will be saved",
    )
    parser.add_argument(
        "--min_class_level",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="Minimum class_level to visualize (inclusive)",
    )
    parser.add_argument(
        "--max_class_level",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="Maximum class_level to visualize (inclusive)",
    )
    parser.add_argument(
        "--ligand2tid_tsv_file",
        type=str,
        required=False,
        default=None,
        help="(Optional) TSV file containing mapping from ligands (molregno) to targets (tid). If given will also gather ligand similarity values per TID.",
    )
    parser.add_argument(
        "--activities_tsv_file",
        type=str,
        required=False,
        default=None,
        help="(Optional) Input TSV file containing active compound information (ligand <-> assay_id mapping). If given will also gather ligand similarity values per ASSAY_ID.",
    )
    args = parser.parse_args()
    return args


def get_overall_dist_stats(
    sim_matrix: np.ndarray,
    similarity_threshold: float = 0.85,
) -> dict:
    total_pairs = len(sim_matrix)
    high_sim_total = (sim_matrix > similarity_threshold).sum()
    results = {"total_pairs": total_pairs, "high_sim_total": high_sim_total}
    return results


def calculate_cluster_probability_for_single_cluster(
    cluster_ligands: set,
    overall_dist_stats: dict,
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    similarity_threshold: float = 0.85,
) -> dict:
    """
    Calculate P(both in cluster | similarity > threshold) for a specific cluster.

    Args:
        cluster_ligands: Set of ligand IDs that belong to this cluster
        all_ligand_ids: List of all ligand IDs to consider
        sim_matrix: 1D array containing similarity values (lower triangular matrix)
        id2idx_map: Mapping from ligand IDs to matrix indices
        similarity_threshold: Threshold for high similarity (default 0.85)

    Returns:
        Dictionary with probability statistics for this cluster
    """
    # Counters for conditional probability calculation
    high_sim_both_in_cluster = 0  # P(both in cluster AND sim > threshold)

    # Additional counters for comprehensive analysis
    total_pairs = overall_dist_stats["total_pairs"]
    high_sim_total = overall_dist_stats["high_sim_total"]
    both_in_cluster_total = 0

    # Iterate through all pairs of ligands
    for id1, id2 in combinations(cluster_ligands, 2):
        # Get similarity value
        i1, i2 = id2idx_map[id1], id2idx_map[id2]
        row = max(i1, i2)
        col = min(i1, i2)
        sim_idx = get_ltm_idx(row, col)
        sim_value = sim_matrix[sim_idx]

        both_in_cluster_total += 1
        if sim_value > similarity_threshold:
            high_sim_both_in_cluster += 1

    # Calculate probabilities
    if high_sim_total == 0:
        conditional_prob = 0.0
    else:
        conditional_prob = high_sim_both_in_cluster / high_sim_total

    # Additional useful statistics
    baseline_prob = both_in_cluster_total / total_pairs if total_pairs > 0 else 0

    results = {
        "conditional_probability": conditional_prob,
        "baseline_probability": baseline_prob,
        "high_sim_both_in_cluster_count": high_sim_both_in_cluster,
        "high_sim_total_count": high_sim_total,
        "total_pairs": total_pairs,
        "both_in_cluster_total": both_in_cluster_total,
        "similarity_threshold": similarity_threshold,
        "enrichment_factor": (
            conditional_prob / baseline_prob if baseline_prob > 0 else float("inf")
        ),
        "cluster_size": len(cluster_ligands),
    }

    return results


def calculate_cluster_probabilities_all_clusters(
    l2c_df: pl.DataFrame,
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    similarity_threshold: float = 0.85,
    group_col: str = "cluster",
    id_col: str = "molregno",
) -> dict:
    """
    Calculate P(both in cluster | similarity > threshold) for each cluster individually.

    Args:
        l2c_df: DataFrame with ligand-to-cluster mappings
        sim_matrix: 1D array containing similarity values (lower triangular matrix)
        id2idx_map: Mapping from ligand IDs to matrix indices
        similarity_threshold: Threshold for high similarity (default 0.85)
        group_col: Column name for cluster assignments
        id_col: Column name for ligand IDs

    Returns:
        Dictionary mapping cluster_id -> probability statistics
    """
    overall_dist_stats = get_overall_dist_stats(sim_matrix, similarity_threshold)
    # Get unique clusters
    unique_clusters = l2c_df[group_col].unique().to_list()

    cluster_results = {}

    for cluster_id in unique_clusters:
        # Get ligands in this cluster
        cluster_ligands = set(
            l2c_df.filter(pl.col(group_col) == cluster_id)[id_col].to_list()
        )

        # Skip clusters with too few ligands (need at least 2 for pairs)
        if len(cluster_ligands) < 2:
            continue

        # Calculate probability statistics for this cluster
        cluster_stats = calculate_cluster_probability_for_single_cluster(
            cluster_ligands,
            overall_dist_stats,
            sim_matrix,
            id2idx_map,
            similarity_threshold,
        )

        cluster_results[cluster_id] = cluster_stats

    return cluster_results


# TODO: use np.linspace to calculate threshold values (take in from args)
def analyze_probability_vs_threshold_per_cluster(
    l2c_df: pl.DataFrame,
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    thresholds: list[float] = [0.75, 0.85, 0.95],
    group_col: str = "cluster",
    id_col: str = "molregno",
) -> pl.DataFrame:
    """
    Analyze how conditional probability changes with threshold for each cluster.

    Returns:
        DataFrame with threshold analysis results per cluster
    """
    results = []

    for threshold in tqdm(thresholds, file=sys.stdout):
        cluster_results = calculate_cluster_probabilities_all_clusters(
            l2c_df, sim_matrix, id2idx_map, threshold, group_col, id_col
        )

        for cluster_id, stats in cluster_results.items():
            result_row = {
                "cluster_id": cluster_id,
                "threshold": threshold,
                "conditional_probability": stats["conditional_probability"],
                "baseline_probability": stats["baseline_probability"],
                "high_sim_pairs": stats["high_sim_both_in_cluster_count"],
                "enrichment_factor": stats["enrichment_factor"],
                "cluster_size": stats["cluster_size"],
            }
            results.append(result_row)

    return pl.DataFrame(results)


def run_probability_analysis(
    l2c_df: pl.DataFrame,
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    tsv_save_path: str,
    group_col: str = "cluster",
    id_col: str = "molregno",
) -> None:
    """
    Run comprehensive probability analysis per cluster and save results.
    """

    # Threshold sweep analysis
    threshold_analysis = analyze_probability_vs_threshold_per_cluster(
        l2c_df, sim_matrix, id2idx_map, group_col=group_col, id_col=id_col
    )

    # Save results
    threshold_analysis.write_csv(tsv_save_path, separator="\t")

    print(f"Saved threshold analysis to: {tsv_save_path}")


def main():
    args = parse_args()
    id_list = load_from_pkl(args.similarity_id_pkl_file)
    N = len(id_list)
    id2idx_map = get_id2idx_map(id_list)
    sim_matrix = np.load(args.similarity_npy_file)
    check_size(sim_matrix, N)

    prob_analysis_dir = args.prob_analysis_dir
    os.makedirs(prob_analysis_dir, exist_ok=True)

    # gather compounds in each cluster and their similarities
    # then save result to pkl file
    for cl in range(args.min_class_level, args.max_class_level + 1):
        l2c_path = get_ligand2cluster_fpath(args.cluster_dir, cl)
        l2c_df = pl.read_csv(l2c_path, separator="\t", columns=["molregno", "cluster"])

        save_path = os.path.join(
            prob_analysis_dir,
            f"per_cluster_threshold_analysis_class_level={cl}.tsv",
        )
        run_probability_analysis(l2c_df, sim_matrix, id2idx_map, save_path)

    # save per tid sim values (if requested)
    if args.ligand2tid_tsv_file is not None:
        ligand2tid_df = pl.read_csv(
            args.ligand2tid_tsv_file, separator="\t", columns=["molregno", "tid"]
        )
        save_path = os.path.join(
            prob_analysis_dir,
            f"per_cluster_threshold_analysis_per_target.tsv",
        )
        run_probability_analysis(
            ligand2tid_df,
            sim_matrix,
            id2idx_map,
            save_path,
            group_col="tid",
        )

    if args.activities_tsv_file is not None:
        assay2tid_df = pl.read_csv(
            args.activities_tsv_file, separator="\t", columns=["molregno", "assay_id"]
        )
        save_path = os.path.join(
            prob_analysis_dir,
            f"per_cluster_threshold_analysis_per_assay.tsv",
        )
        run_probability_analysis(
            assay2tid_df,
            sim_matrix,
            id2idx_map,
            save_path,
            group_col="assay_id",
        )


if __name__ == "__main__":
    main()
