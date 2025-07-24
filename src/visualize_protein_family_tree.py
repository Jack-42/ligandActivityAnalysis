"""
@author Jack Ringer
Date: 6/23/2025
Description:
Code to visualize a phylogenic tree of the chosen
protein family (e.g., the kinome).
Currently this generates just the family tree image without any legend or title - these have to be added separately.
I'll admit this script is a mess and could likely be improved - was just generated quick and dirty w/ Claude (mostly).
"""

import argparse

import polars as pl
import seaborn as sns
from bigtree import polars_to_tree_by_relation, tree_to_newick
from ete4 import Tree
from ete4.treeview import NodeStyle, TreeStyle


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
    parser.add_argument(
        "--major_family_short_names",
        type=lambda s: [item.upper() for item in s.split(",")],
        required=True,
        help="List of short_name denoting major families (ignores case). Will color the tree by these provided families with the given pallette. See logs for how each color maps to each major family (trying to add legends with ETE4 is too much of a pain).",
    )
    parser.add_argument(
        "--pallette",
        type=str,
        default="hls",
        help="Seaborn pallette to use for coloring tree by major family (sns.color_palette(pallette)). Will print out hex colors from this pallette to support manually adding in legend.",
    )
    args = parser.parse_args()
    return args


def get_major_family_colors(major_family_short_names: list[str], pallette: str):
    return dict(
        zip(
            major_family_short_names,
            sns.color_palette(
                pallette, n_colors=len(major_family_short_names)
            ).as_hex(),
        )
    )


def find_major_family_for_node(node, df, major_families):
    # Create a mapping from protein_class_id to parent_id and short_name
    id_to_parent = dict(zip(df["protein_class_id"], df["parent_id"]))
    id_to_short_name = dict(zip(df["protein_class_id"], df["short_name"]))

    # Extract protein_class_id from node name (assuming it's in the format like "119_Pdhk_5")
    try:
        node_id = int(node.name.split("_")[0])
    except (ValueError, IndexError):
        return None

    current_id = node_id
    visited = set()

    while current_id and current_id not in visited:
        visited.add(current_id)
        current_short_name = id_to_short_name.get(current_id)

        if current_short_name in major_families:
            return current_short_name

        # Move to parent
        current_id = id_to_parent.get(current_id)

    return None


def customize_tree_appearance(tree: Tree, df: pl.DataFrame, major_family_colors: dict):
    """Customize the appearance of the tree with colors and labels"""
    major_families = set(major_family_colors.keys())

    for node in tree.traverse():
        # Create node style
        nstyle = NodeStyle()

        major_family = find_major_family_for_node(node, df, major_families)

        group_color = "#CCCCCC"
        if major_family:
            group_color = major_family_colors[major_family]

        nstyle["hz_line_color"] = group_color
        nstyle["vt_line_color"] = group_color
        nstyle["hz_line_width"] = 1
        nstyle["vt_line_width"] = 1
        nstyle["size"] = 0

        # Set node style
        node.set_style(nstyle)


def main():
    args = parse_args()
    df = pl.read_csv(args.family_details_tsv_file, separator="\t")
    df = df.with_columns((pl.col("short_name").str.to_uppercase()))

    # Build tree
    t = polars_to_tree_by_relation(
        df, child_col="protein_class_id", parent_col="parent_id"
    )
    nwt = tree_to_newick(t, attr_list=["short_name", "class_level"])
    nwt = nwt + ";"
    et = Tree(nwt, parser=1)

    # Customize tree appearance
    major_family_colors = get_major_family_colors(
        args.major_family_short_names, args.pallette
    )

    customize_tree_appearance(et, df.to_pandas(), major_family_colors)

    # Set up circular tree style
    circular_style = TreeStyle()
    circular_style.mode = "c"  # circular mode
    circular_style.scale = 20
    circular_style.show_leaf_name = False
    circular_style.show_branch_length = False
    circular_style.show_branch_support = False
    circular_style.show_scale = False

    # print legend
    print("LEGEND (short_name to hex color):", major_family_colors)

    # Render the tree
    et.render(
        args.out_path, w=1024, h=1024, units="px", tree_style=circular_style, dpi=300
    )


if __name__ == "__main__":
    main()
