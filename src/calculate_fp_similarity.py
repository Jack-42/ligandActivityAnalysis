"""
@author Jack Ringer
Date: 7/3/2025
Description:
Calculate similarity matrix from the output
of generate_fingerprints.py (Python dictionary).
Currently calculates Tanimoto coeffs, but could be
adapted to other similarity metrics in the future.
"""

import argparse

import numpy as np
from rdkit.DataManip.Metric.rdMetricMatrixCalc import GetTanimotoSimMat

from utils.io import load_from_pkl, save_to_pkl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate lower-triangular similarity matrix from input fingerprints and save to .npy file",
        epilog="",
    )
    parser.add_argument(
        "--fingerprints_pkl_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input Pickle file containing compound ID mapped to fingerprint",
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
        "--max_n_compounds",
        type=int,
        default=20_000,
        help="Maximum number of fingerprints to be processed. Included as a failsafe to help warn users of extreme computation/storage cost",
    )
    args = parser.parse_args()
    return args


def assertion_msg(d: dict, max_n_compounds: int) -> str:
    return f"""
    Given number of compounds was {len(d)} when max_n_compounds={max_n_compounds}. 
    Note that creating a full similarity matrix is O(N^2) in terms of both computational cost and storage, 
    so please make sure to set a reasonable limit on max_n_compounds. 
    If you have a very large number of compounds then you may want to use other methods 
    (e.g., sampling random pairs of molecules).
    """


def main():
    args = parse_args()
    assert args.similarity_npy_file.endswith(".npy")
    d = load_from_pkl(args.fingerprints_pkl_file)
    assert len(d) < args.max_n_compounds, assertion_msg(d, args.max_n_compounds)
    ids, fingerprints = zip(*d.items())
    sim_matrix = GetTanimotoSimMat(fingerprints)
    np.save(args.similarity_npy_file, sim_matrix)
    save_to_pkl(ids, args.similarity_id_pkl_file)


if __name__ == "__main__":
    main()
