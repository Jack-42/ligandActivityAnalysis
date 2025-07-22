"""
@author Jack Ringer
Date: 7/22/2025
Description:
This script will compare the overall distribution of similarity values (- the group)
to individual group distributions using the Mann-Whitney U-test.
Note that I've written this script to deal with the distributions
observed for Kinase data. For other protein families a different
statistical test/approach may be more appropriate.
"""

import argparse
from math import comb

import numpy as np
import polars as pl
from scipy.stats import mannwhitneyu

from utils.constants import get_cluster2sim_fpath, get_ligand2cluster_fpath
from utils.io import load_from_pkl
from utils.ltm import check_size, get_all_ltm_indices_from_idx, get_id2idx_map


def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform Mann-Whitney U-test comparing overall distribution of similarity values to subgroups/families",
        epilog="",
    )
    parser.add_argument(
        "--similarity_npy_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output .npy file containing lower-triangular similarity matrix",
    )
    parser.add_argument(
        "--cluster_dir",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input directory where clustered ligands and targets were saved to by assign_family_clusters.py",
    )
    parser.add_argument(
        "--cluster2sim_dir",
        type=str,
        required=True,
        help="Input directory where .pkl files containing similarity values gathered per class_level + cluster",
    )
    parser.add_argument(
        "--similarity_id_pkl_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output .pkl file containing compound IDs (in order used by similarity_npy_file)",
    )
    parser.add_argument(
        "--family_details_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing protein_class_id and additional information (pref_name, short_name, class_level)",
    )
    parser.add_argument(
        "--class_level",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="The class_level to run the test on",
    )
    parser.add_argument(
        "--mann_whitney_utest_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing results from running the mann-whitney u-test on each cluster",
    )
    args = parser.parse_args()
    return args


def sanity_checks(ligand_cluster_df: pl.DataFrame, cluster_dict: dict):
    assert len(cluster_dict.keys()) == ligand_cluster_df["cluster"].n_unique()
    for cluster_tup, cluster_df in ligand_cluster_df.group_by("cluster"):
        cluster_id = cluster_tup[0]
        dict_n_pairs = len(cluster_dict[cluster_id])
        df_n_pairs = comb(cluster_df["molregno"].n_unique(), 2)
        assert (
            dict_n_pairs == df_n_pairs
        ), f"Mismatch with cluster={cluster_id} ({dict_n_pairs} != {df_n_pairs})"


def main():
    args = parse_args()
    id_list = load_from_pkl(args.similarity_id_pkl_file)
    N = len(id_list)
    id2idx_map = get_id2idx_map(id_list)
    sim_matrix = np.load(args.similarity_npy_file)
    check_size(sim_matrix, N)

    family_info_df = pl.read_csv(args.family_details_tsv_file, separator="\t")

    cluster_dict = load_from_pkl(
        get_cluster2sim_fpath(args.cluster2sim_dir, args.class_level)
    )

    ligand_cluster_df = pl.read_csv(
        get_ligand2cluster_fpath(args.cluster_dir, args.class_level), separator="\t"
    )

    ligand_cluster_df = ligand_cluster_df.join(family_info_df, on="protein_class_id")
    sanity_checks(ligand_cluster_df, cluster_dict)

    cluster2result_dict = {
        "cluster": [],
        "cluster_size": [],
        "other_size": [],
        "cluster_median": [],
        "other_median": [],
        "U1": [],
        "U2": [],
        "p_val": [],
        "corrected_p_val": [],  # using a Bonferroni correction to avoid multiple comparisons problem
    }
    n_clusters = len(cluster_dict)
    for cluster in cluster_dict:
        cluster_sim_values = cluster_dict[cluster]
        # need to remove similarity values computed wrt to all ligands belonging to this cluster for test to be valid
        # note we can't just add together other clusters to get other_sim_values due to multitarget/group activity of some ligands
        cluster_ligands = (
            ligand_cluster_df.filter(ligand_cluster_df["cluster"] == cluster)[
                "molregno"
            ]
            .unique()
            .to_list()
        )
        cluster_ligand_indices = [id2idx_map[lig_id] for lig_id in cluster_ligands]
        cluster_indices = []
        for i in cluster_ligand_indices:
            ligand_indices = get_all_ltm_indices_from_idx(i, N)
            assert len(ligand_indices) == (
                N - 1
            ), f"Each ligand would have been compared against N-1 = {N-1} other ligands, but found len(ligand_indices)={len(ligand_indices)}"
            cluster_indices += ligand_indices

        other_sim_values = np.delete(sim_matrix, cluster_indices)

        # "greater" because we want to know if the distribution of cluster_sim_values is stochastically greater than other_sim_values
        U1, p_val = mannwhitneyu(
            cluster_sim_values, other_sim_values, alternative="greater"
        )
        # see notes of: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html
        cluster_size = len(cluster_sim_values)
        other_size = len(other_sim_values)
        U2 = cluster_size * other_size - U1
        cluster2result_dict["cluster"] += [cluster]
        cluster2result_dict["cluster_size"] += [cluster_size]
        cluster2result_dict["other_size"] += [other_size]
        cluster2result_dict["cluster_median"] += [np.median(cluster_sim_values)]
        cluster2result_dict["other_median"] += [np.median(other_sim_values)]
        cluster2result_dict["U1"] += [U1]
        cluster2result_dict["U2"] += [U2]
        cluster2result_dict["p_val"] += [p_val]
        cluster2result_dict["corrected_p_val"] += [min(p_val * n_clusters, 1.0)]

    cluster2name_df = (
        ligand_cluster_df[["cluster", "short_name"]].unique().sort(by="cluster")
    )
    result_df = pl.from_dict(cluster2result_dict)
    result_df = result_df.join(cluster2name_df, on="cluster")
    result_df = result_df[
        [
            "cluster",
            "short_name",
            "cluster_size",
            "other_size",
            "cluster_median",
            "other_median",
            "U1",
            "U2",
            "p_val",
            "corrected_p_val",
        ]
    ]
    result_df = result_df.sort(by="cluster")
    result_df.write_csv(args.mann_whitney_utest_tsv_file, separator="\t")


if __name__ == "__main__":
    main()
