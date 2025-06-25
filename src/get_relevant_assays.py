"""
@author Jack Ringer
Date: 6/16/2025
Description:
From a given TSV file containing target IDs (TID), fetch
all assays exploring activity of the given targets.
"""

import argparse
from typing import Optional

import polars as pl
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
    parser.add_argument(
        "--doc_type",
        type=str,
        required=False,
        default=None,
        help="(Optional) Filter to only include assays associated with the given doc_type",
    )
    parser.add_argument(
        "--confidence_score",
        type=int,
        required=False,
        default=None,
        help="(Optional) Filter by given confidence_score",
    )
    parser.add_argument(
        "--exclude_variants",
        action=argparse.BooleanOptionalAction,
        help="Exclude assays which target a variant. Based on criteria described by: https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049",
    )
    args = parser.parse_args()
    return args


def get_assays(
    cursor: DictCursor,
    tids: list[int],
    assay_type: Optional[str] = None,
    confidence_score: Optional[int] = None,
    doc_type: Optional[str] = None,
    exclude_variants: bool = False,
) -> list[int]:
    query = sql.SQL(
        """
        SELECT assay_id, tid
        FROM assays ass
        JOIN docs d ON ass.doc_id = d.doc_id
        WHERE tid IN {tids}
        AND ({assay_type} IS NULL OR assay_type = {assay_type})
        AND ({confidence_score} IS NULL OR confidence_score = {confidence_score})
        AND ({doc_type} IS NULL OR d.doc_type = {doc_type})
        AND (NOT {exclude_variants} OR (variant_id is NULL AND NOT(LOWER(description) LIKE '%mutant%' OR LOWER(description) LIKE '%mutation%' OR LOWER(description) LIKE '%variant%')));
        """
    ).format(
        tids=sql.Literal(tuple(tids)),
        assay_type=sql.Literal(assay_type),
        confidence_score=sql.Literal(confidence_score),
        doc_type=sql.Literal(doc_type),
        exclude_variants=sql.Literal(exclude_variants),
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

    target_df = pl.read_csv(args.target_tsv_file, separator="\t")
    tids = list(set(target_df["tid"]))
    assays = get_assays(
        cursor,
        tids,
        args.assay_type,
        args.confidence_score,
        args.doc_type,
        args.exclude_variants,
    )
    write_to_tsv(args.assay_tsv_file, assays, ["assay_id", "tid"])

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
