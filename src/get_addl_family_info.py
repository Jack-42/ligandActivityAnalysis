"""
@author Jack Ringer
Date: 6/18/2025
Description:
Fetch additional information re protein families from the
ChEMBL DB.
Right now this logic could just be part of get_protein_targets.py, but keeping it
separate allows for greater flexibility in future (if need be).
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
        description="Get additional information on the given protein families (protein_class_id).",
        epilog="",
    )
    add_chembl_db_args(parser)
    parser.add_argument(
        "--family_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing family relationship information from protein_classification table",
    )
    parser.add_argument(
        "--family_details_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing protein_class_id and additional information (pref_name, short_name, class_level)",
    )
    args = parser.parse_args()
    return args


def get_addl_info(cursor: DictCursor, protein_class_ids: list[int]):
    query = sql.SQL(
        """
        SELECT protein_class_id,parent_id,pref_name,short_name,class_level
        FROM protein_classification
        WHERE protein_class_id IN {protein_class_ids};
        """
    ).format(protein_class_ids=sql.Literal(tuple(protein_class_ids)))

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

    df = pl.read_csv(args.family_tsv_file, separator="\t")
    protein_class_ids = list(df["protein_class_id"])
    family_details = get_addl_info(cursor, protein_class_ids)
    write_to_tsv(
        args.family_details_tsv_file,
        family_details,
        ["protein_class_id", "parent_id", "pref_name", "short_name", "class_level"],
    )

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
