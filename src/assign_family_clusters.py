"""
@author Jack Ringer
Date: 6/30/2025
Description:
Group together ligands which are active against the same
family / protein.
"""

import argparse

import polars as pl
from bigtree import find_path, polars_to_tree_by_relation

from utils.constants import get_ligand2cluster_fpath, get_target2cluster_fpath


def parse_args():
    parser = argparse.ArgumentParser(
        description="Assign cluster(s) to a given ligands. Will create separate files for each value of class_level",
        epilog="",
    )
    parser.add_argument(
        "--activities_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing active compound information",
    )
    parser.add_argument(
        "--assay_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing TID and assay_id",
    )
    parser.add_argument(
        "--target_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing TID and protein_class_id for each protein",
    )
    parser.add_argument(
        "--family_details_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing protein_class_id and additional information (pref_name, short_name, class_level)",
    )
    parser.add_argument(
        "--cluster_dir",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output directory to save clustered ligands and targets to",
    )
    parser.add_argument(
        "--ligand2tid_tsv_file",
        type=str,
        required=False,
        default=None,
        help="(Optional) Output TSV file containing mapping from ligands (molregno) to targets (tid)",
    )
    parser.add_argument(
        "--min_class_level",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="Minimum class_level to cluster by (inclusive)",
    )
    parser.add_argument(
        "--max_class_level",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="Maximum class_level to cluster by (inclusive)",
    )
    args = parser.parse_args()
    return args


def get_protein_class_id_of_class_level(leaf_node, class_level: int) -> int:
    if leaf_node.class_level < class_level:
        return None
    node = leaf_node
    while node.class_level > class_level:
        node = node.parent
    return node.name  # will be protein_class_id based on main() construction


def main():
    args = parse_args()
    target_df = pl.read_csv(
        args.target_tsv_file, separator="\t", columns=["tid", "protein_class_id"]
    )
    fam_df = pl.read_csv(
        args.family_details_tsv_file,
        separator="\t",
        columns=["protein_class_id", "parent_id", "class_level"],
    )
    activity_df = pl.read_csv(
        args.activities_tsv_file, separator="\t", columns=["molregno", "assay_id"]
    )
    assay_df = pl.read_csv(
        args.assay_tsv_file, separator="\t", columns=["assay_id", "tid"]
    )

    # map active ligands to targets
    joined = activity_df.join(assay_df, on="assay_id", how="inner")
    mol2tid_df = joined.select(["molregno", "tid"]).unique()
    assert len(mol2tid_df["molregno"].unique()) == len(activity_df["molregno"].unique())

    # load family relationships
    tree = polars_to_tree_by_relation(
        fam_df, child_col="protein_class_id", parent_col="parent_id"
    )
    target_df = target_df.with_columns(
        target_df["protein_class_id"]
        .map_elements(lambda pid: find_path(tree, f"/{pid}"), return_dtype=pl.Object)
        .alias("fam_path")
    )
    assert not (target_df["fam_path"].is_null().any())

    # cluster ligands by protein family info
    # NOTE: a single ligand can belong to multiple clusters because of:
    # 1) ligand being active against > 1 target
    # 2) ligand active against a single target, but that target has multiple classifications
    # map ligands to active target(s) + family info
    for cl in range(args.min_class_level, args.max_class_level + 1):
        tid2cluster_df = target_df.with_columns(
            target_df["fam_path"]
            .map_elements(
                lambda n: get_protein_class_id_of_class_level(n, cl), return_dtype=int
            )
            .alias("protein_class_id")
        )
        # null indicates no classification at this class_level
        tid2cluster_df = tid2cluster_df.drop_nulls(subset=["protein_class_id"])
        # help make it explicit that we're clustering by protein_class_id
        unique_clusters = tid2cluster_df["protein_class_id"].unique().sort()
        cluster_mapping = {cluster: i + 1 for i, cluster in enumerate(unique_clusters)}
        tid2cluster_df = tid2cluster_df.with_columns(
            pl.col("protein_class_id").replace(cluster_mapping).alias("cluster")
        )
        tid2cluster_df = tid2cluster_df[["tid", "protein_class_id", "cluster"]]
        mol2cluster_df = mol2tid_df.join(tid2cluster_df, on="tid", how="inner")
        mol2cluster_df = mol2cluster_df.sort(by="molregno")
        tid2cluster_df = tid2cluster_df.sort(by="tid")
        mol_save_path = get_ligand2cluster_fpath(args.cluster_dir, cl)
        tid_save_path = get_target2cluster_fpath(args.cluster_dir, cl)
        mol2cluster_df.write_csv(mol_save_path, separator="\t")
        # saving tid2cluster_df bc there may be some targets which have no active ligands
        tid2cluster_df.write_csv(tid_save_path, separator="\t")

    # save mol2tid dataframe
    if args.ligand2tid_tsv_file is not None:
        mol2tid_df = mol2tid_df.sort(by="molregno")
        mol2tid_df.write_csv(args.ligand2tid_tsv_file, separator="\t")


if __name__ == "__main__":
    main()
