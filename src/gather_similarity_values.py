"""
@author Jack Ringer
Date: 7/7/2025
Description:
Visualize distributions of similarity values between multiple
clusters.
"""

import argparse
from itertools import combinations

import numpy as np
import polars as pl
from rdkit import DataStructs  # for validation

from utils.constants import (
    get_assay2sim_fpath,
    get_cluster2sim_fpath,
    get_ligand2cluster_fpath,
    get_tid2sim_fpath,
)
from utils.io import load_from_pkl, save_to_pkl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Gather similarity values contained in lower triangular matrix (LTM) based on given cluster files",
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
        "--cluster2sim_dir",
        type=str,
        required=True,
        help="Output directory where .pkl files containing similarity values gathered per class_level + cluster",
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
        "--fingerprints_pkl_file",
        type=str,
        default=None,
        help="(Optional) Input Pickle file containing compound ID mapped to fingerprint. If given will be used to validate the logic of this script.",
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


def get_ltm_idx(row: int, col: int) -> int:
    """
    In the LTM (e.g., from GetTanimotoSimMat), find the
    coefficient value between the given row/column.

    :param int row: index of molecule 1 in id_list ("row" in matrix)
    :param int col: index of molecule 2 in id_list ("column" in matrix)
    :return int: index in the LTM with coeff value
    """
    assert col < row
    offset = (row * (row - 1)) // 2
    return offset + col


def run_check(fetched_sim: float, fp_dict: dict, id1: int, id2: int, i1: int, i2: int):
    true_sim = DataStructs.TanimotoSimilarity(fp_dict[id1], fp_dict[id2])
    assert np.isclose(
        fetched_sim, true_sim
    ), f"ids={(id1, id2)}\nindices={(i1, i2)}\nfetched_sim={fetched_sim}\ntrue_sim={true_sim}\nLogic for get_ltm_idx is incorrect, check logic."


def get_similarity_values(
    ligand_ids: list[int],
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    validate: bool,
    fp_dict: dict = None,
) -> list[float]:
    sim_values = []
    for id1, id2 in combinations(ligand_ids, 2):
        i1, i2 = id2idx_map[id1], id2idx_map[id2]
        row = max(i1, i2)
        col = min(i1, i2)
        sim_idx = get_ltm_idx(row, col)
        sim_value = sim_matrix[sim_idx]
        sim_values.append(sim_value)
        if validate:
            run_check(sim_value, fp_dict, id1, id2, i1, i2)
    return sim_values


def gather_and_save_cluster_sim_values(
    l2c_df: pl.DataFrame,
    sim_matrix: np.ndarray,
    id2idx_map: dict,
    save_path: str,
    validate: bool = False,
    fp_dict: dict = None,
    group_col: str = "cluster",
    id_col: str = "molregno",
) -> None:
    cluster2simvalues = {}
    for cluster_tup, cluster_df in l2c_df.group_by(group_col):
        cluster_ligand_ids = cluster_df[id_col].unique().to_list()
        sim_values = get_similarity_values(
            cluster_ligand_ids, sim_matrix, id2idx_map, validate, fp_dict
        )
        cluster_id = cluster_tup[0]  # polars quirk
        cluster2simvalues[cluster_id] = sim_values
    save_to_pkl(cluster2simvalues, save_path)


def main():
    args = parse_args()
    id_list = load_from_pkl(args.similarity_id_pkl_file)
    N = len(id_list)
    id2idx_map = {cpd_id: i for i, cpd_id in enumerate(id_list)}
    sim_matrix = np.load(args.similarity_npy_file)
    assert (
        len(sim_matrix.shape) == 1
    ), f"Expected LTM to be a 1D numpy array, given shape: {sim_matrix.shape}"
    assert len(sim_matrix) == (
        N * (N - 1) // 2
    ), "Mismatch between length of given id_list and LTM, check input"

    # check if we do extra validation
    fp_dict = {}
    validate = False
    if args.fingerprints_pkl_file is not None:
        fp_dict = load_from_pkl(args.fingerprints_pkl_file)
        validate = True

    # gather compounds in each cluster and their similarities
    # then save result to pkl file
    for cl in range(args.min_class_level, args.max_class_level + 1):
        l2c_path = get_ligand2cluster_fpath(args.cluster_dir, cl)
        l2c_df = pl.read_csv(l2c_path, separator="\t", columns=["molregno", "cluster"])
        save_path = get_cluster2sim_fpath(args.cluster2sim_dir, cl)
        gather_and_save_cluster_sim_values(
            l2c_df, sim_matrix, id2idx_map, save_path, validate, fp_dict
        )

    # save per tid sim values (if requested)
    if args.ligand2tid_tsv_file is not None:
        ligand2tid_df = pl.read_csv(
            args.ligand2tid_tsv_file, separator="\t", columns=["molregno", "tid"]
        )
        save_path = get_tid2sim_fpath(args.cluster2sim_dir)
        gather_and_save_cluster_sim_values(
            ligand2tid_df,
            sim_matrix,
            id2idx_map,
            save_path,
            validate,
            fp_dict,
            group_col="tid",
        )

    if args.activities_tsv_file is not None:
        assay2tid_df = pl.read_csv(
            args.activities_tsv_file, separator="\t", columns=["molregno", "assay_id"]
        )
        save_path = get_assay2sim_fpath(args.cluster2sim_dir)
        gather_and_save_cluster_sim_values(
            assay2tid_df,
            sim_matrix,
            id2idx_map,
            save_path,
            validate,
            fp_dict,
            group_col="assay_id",
        )


if __name__ == "__main__":
    main()
