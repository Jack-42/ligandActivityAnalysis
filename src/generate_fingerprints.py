"""
@author Jack Ringer
Date: 6/26/2025
Description: Generate molecular fingerprints from
input file containing mol SMILES.
"""

import argparse

import polars as pl
from rdkit import Chem
from rdkit.Chem import AllChem


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate molecular features from SMILES tsv file and save the features to a Parquet file.",
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
        "--fingerprints_parquet_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output Parquet file containing compound ID mapped to fingerprint",
    )
    args = parser.parse_args()
    return args


def smiles_to_fp(smiles: str, radius: int = 2, nBits: int = 2048) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
    return fp.ToBitString()


def main():
    args = parse_args()
    df = pl.read_csv(args.compound_structures_tsv_file, separator="\t")
    df = df.with_columns(
        pl.col("canonical_smiles").map_elements(smiles_to_fp).alias("fingerprint")
    )
    assert not (
        df["fingerprint"].is_null().any()
    ), "Null fingerprint was generated, check input TSV file and ensure all SMILEs are valid"
    df[["molregno", "fingerprint"]].write_parquet(args.fingerprints_parquet_file)


if __name__ == "__main__":
    main()
