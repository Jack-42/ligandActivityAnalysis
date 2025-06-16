"""
@author Jack Ringer
Date: 6/16/2025
Description:
Simple script to plot a histogram of the number of unique
assays per target.
"""

import argparse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(
        description="Return the TargetIDs (TIDs) of proteins belonging to the given family, as well as protein family relation information (e.g., if the given family has subfamilies).",
        epilog="",
    )
    parser.add_argument(
        "--assay_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing TID and assay_id",
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


def main():
    args = parse_args()
    df = pd.read_csv(args.assay_tsv_file, sep="\t")
    n_a_per_t = df["tid"].value_counts().values
    sns.histplot(n_a_per_t)
    plt.title("Number of Unique Assay Records per Target")
    plt.xlabel("N. Unique Assays")
    plt.ylabel("Frequency")
    plt.savefig(args.out_path, dpi=300)


if __name__ == "__main__":
    main()
