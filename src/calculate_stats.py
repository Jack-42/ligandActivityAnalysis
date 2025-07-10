"""
@author Jack Ringer
Date: 7/10/2025
Description:
Compute statistics from the similarity values of the overall
population (e.g., all protein kinase ligands) and the
sub-clusters (e.g., kinase groups).
"""

import argparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate statistics (mean, median, std) from the provided files",
        epilog="",
    )
    args = parser.parse_args()
    return args
