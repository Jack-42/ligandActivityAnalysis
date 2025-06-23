"""
@author Jack Ringer
Date: 6/23/2025
Description:
Code to visualize a phylogenic tree of the chosen
protein family (e.g., the kinome).
"""

import argparse

import polars as pl
from bigtree import polars_to_tree_by_relation, tree_to_newick
from ete4 import Tree
from ete4.treeview import TreeStyle


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create visualization of a protein family's phylogenic tree based on the provided TSV file.",
        epilog="",
    )
    parser.add_argument(
        "--family_details_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing protein_class_id and additional information (pref_name, short_name, class_level)",
    )
    parser.add_argument(
        "--out_path",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output path to save png to",
    )
    args = parser.parse_args()
    return args


# TODO: add addl tree customization, color and show regions of important subfamilies
def main():
    args = parse_args()
    df = pl.read_csv(args.family_details_tsv_file, separator="\t")
    t = polars_to_tree_by_relation(
        df, child_col="protein_class_id", parent_col="parent_id"
    )
    nwt = tree_to_newick(t, attr_list=["short_name", "class_level"])
    nwt = nwt + ";"
    et = Tree(nwt, parser=1)
    circular_style = TreeStyle()
    circular_style.mode = "c"  # draw tree in circular mode
    circular_style.scale = 20
    et.render(args.out_path, w=1024, units="mm", tree_style=circular_style)


if __name__ == "__main__":
    main()
