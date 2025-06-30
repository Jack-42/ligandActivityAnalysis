"""
@author Jack Ringer
Date: 6/26/2025
Description: Generate molecular fingerprints from
input file containing mol SMILES.
"""

import argparse

import numpy as np
import polars as pl
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFingerprintGenerator

from utils.io import save_to_pkl


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate molecular features from SMILES tsv file and save the features to a Pickle file.",
        epilog="",
    )
    parser.add_argument(
        "--compound_structures_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing compound ID and SMILES",
    )
    parser.add_argument(
        "--fingerprints_pkl_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output Pickle file containing compound ID mapped to fingerprint",
    )
    args = parser.parse_args()
    return args


def smiles_to_fp(smiles: str, fp_gen) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = fp_gen.GetFingerprint(mol)
    return fp


def main():
    args = parse_args()
    df = pl.read_csv(args.compound_structures_tsv_file, separator="\t")
    mfpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    df = df.with_columns(
        pl.col("canonical_smiles")
        .map_elements(lambda x: smiles_to_fp(x, mfpgen))
        .alias("fingerprint")
    )
    assert not (
        df["fingerprint"].is_null().any()
    ), "Null fingerprint was generated, check input TSV file and ensure all SMILEs are valid"
    fp_dict = dict(zip(df["molregno"], df["fingerprint"]))
    save_to_pkl(fp_dict, args.fingerprints_pkl_file)


if __name__ == "__main__":
    main()
