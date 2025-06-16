"""
@author Jack Ringer
Date: 6/16/2025
Description:
From a given TSV file containing target IDs (TID), fetch
all assays exploring activity of the given targets.
"""

import argparse
from typing import Optional

import pandas as pd
import psycopg2
from psycopg2 import sql
from psycopg2.extras import DictCursor

from utils.args import add_chembl_db_args
from utils.io import write_to_tsv


def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch assays matching the given targets.",
        epilog="",
    )
    add_chembl_db_args(parser)
    parser.add_argument(
        "--target_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing TID and protein_class_id for each protein",
    )
    parser.add_argument(
        "--assay_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing TID and assay_id",
    )
    parser.add_argument(
        "--assay_type",
        type=str,
        required=False,
        default=None,
        help="(Optional) Filter by given assay_type",
    )
    args = parser.parse_args()
    return args


def get_assays(
    cursor: DictCursor, tids: list[int], assay_type: Optional[str] = None
) -> list[int]:
    query = sql.SQL(
        """
        SELECT assay_id, tid
        FROM assays
        WHERE tid IN {tids}
        AND ({assay_type} IS NULL OR assay_type = {assay_type});
        """
    ).format(
        tids=sql.Literal(tuple(tids)),
        assay_type=sql.Literal(assay_type),
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

    target_df = pd.read_csv(args.target_tsv_file, sep="\t")
    tids = list(set(target_df["tid"]))
    assays = get_assays(cursor, tids, args.assay_type)
    write_to_tsv(args.assay_tsv_file, assays, ["assay_id", "tid"])

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
