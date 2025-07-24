"""
@author Jack Ringer
Date: 6/25/2025
Description:
Get structures (SMILES and Inchi) of given compounds from ChEMBL database
and save to a TSV file.
"""

import argparse

import polars as pl
import psycopg2
from psycopg2 import sql
from psycopg2.extras import DictCursor

from utils.args import add_chembl_db_args
from utils.io import write_to_tsv


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch the IDs of compounds deemed active in the given assays.",
        epilog="",
    )
    add_chembl_db_args(parser)
    parser.add_argument(
        "--activities_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing active compound information",
    )
    parser.add_argument(
        "--compound_structures_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing compound ID, Inchi, and SMILES",
    )
    args = parser.parse_args()
    return args


def get_compound_structures(
    cursor: DictCursor,
    molregno_list: list[int],
) -> list[list]:
    query = sql.SQL(
        """
        SELECT cs.molregno, md.chembl_id, cs.standard_inchi, cs.canonical_smiles
        FROM compound_structures cs
        JOIN molecule_dictionary md on cs.molregno = md.molregno
        WHERE (cs.molregno in {molregno_list});
        """
    ).format(
        molregno_list=sql.Literal(tuple(molregno_list)),
    )
    cursor.execute(query)
    result = cursor.fetchall()
    return result


def main():
    args = parse_args()
    connection = psycopg2.connect(
        dbname=args.db_name,
        host=args.db_host,
        user=args.db_user,
        password=args.db_password,
        port=args.db_port,
        cursor_factory=DictCursor,
    )
    connection.set_session(readonly=True)
    cursor = connection.cursor()

    df = pl.read_csv(args.activities_tsv_file, separator="\t")
    molregno_list = list(set(df["molregno"]))
    structure_info = get_compound_structures(
        cursor,
        molregno_list,
    )
    header = ["molregno", "chembl_id", "standard_inchi", "canonical_smiles"]
    write_to_tsv(args.compound_structures_tsv_file, structure_info, header)

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
