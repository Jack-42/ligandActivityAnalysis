"""
@author Jack Ringer
Date: 6/30/2025
Description:
Group together ligands which are active against the same
family / protein.
"""

import argparse

import polars as pl
from bigtree import polars_to_tree_by_relation


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
        "--output_dir",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output directory to save clustered ligands to",
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


def main():
    args = parse_args()
    activity_df = pl.read_csv(
        args.activities_tsv_file, separator="\t", columns=["molregno", "assay_id"]
    )
    assay_df = pl.read_csv(
        args.assay_tsv_file, separator="\t", columns=["assay_id", "tid"]
    )
    target_df = pl.read_csv(
        args.target_tsv_file, separator="\t", columns=["tid", "protein_class_id"]
    )
    fam_df = pl.read_csv(
        args.family_details_tsv_file,
        separator="\t",
        columns=["protein_class_id", "parent_id" "class_level"],
    )
    t2f_df = target_df.join(fam_df, on="")
    tree = polars_to_tree_by_relation(
        fam_df, child_col="protein_class_id", parent_col="parent_id"
    )
    # TODO: use tree to help find parent / child relationships


if __name__ == "__main__":
    main()
