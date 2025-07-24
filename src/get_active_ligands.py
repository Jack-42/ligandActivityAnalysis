"""
@author Jack Ringer
Date: 6/16/2025
Description:
From a given list of assays, get all ligands which
meet the given criteria for being "active".
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
        description="Fetch the IDs of ligands (i.e., compounds) deemed active in the given assays.",
        epilog="",
    )
    add_chembl_db_args(parser)
    parser.add_argument(
        "--assay_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Input TSV file containing TID and assay_id",
    )
    parser.add_argument(
        "--activities_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing active compound information",
    )
    parser.add_argument(
        "--pchembl_min_value",
        type=float,
        default=5.0,
        help="Minimum activity threshold (defined by pchembl_value)",
    )
    parser.add_argument(
        "--min_mw",
        type=float,
        default=200.0,
        help="Minimum full molecular weight of selected compounds (in Da)",
    )
    parser.add_argument(
        "--max_mw",
        type=float,
        default=900.0,
        help="Maximum full molecular weight of selected compounds (in Da)",
    )
    parser.add_argument(
        "--structural_alert_set_ids",
        type=lambda s: [int(item) for item in s.split(",")],
        default=None,
        help="(Optional) One or more alert_set_ids to filter on (comma-separated)",
    )
    args = parser.parse_args()
    return args


def get_active_mol_info(
    cursor: DictCursor,
    assay_ids: list[int],
    pchembl_min_value: int,
    min_mw: float,
    max_mw: float,
    structural_alert_set_ids: Optional[list[int]] = None,
) -> list[list]:
    query = sql.SQL(
        """
        SELECT md.molregno,md.chembl_id,a.assay_id,a.standard_type,a.pchembl_value
        FROM activities a
        JOIN molecule_dictionary md ON a.molregno = md.molregno
        JOIN compound_properties cp ON a.molregno = cp.molregno
        JOIN compound_structural_alerts csa ON a.molregno = csa.molregno
        JOIN structural_alerts sa ON csa.alert_id = sa.alert_id
        WHERE a.assay_id IN {assay_ids}
        AND a.pchembl_value >= {pchembl_min_value}
        AND a.relation = '='
        AND a.standard_flag = 1
        AND md.structure_type = 'MOL'
        AND (a.data_validity_comment IS NULL OR a.data_validity_comment = 'Manually validated')
        AND (a.potential_duplicate = 0)
        AND (cp.full_mwt >= {min_mw})
        AND (cp.full_mwt <= {max_mw})
        AND (sa.alert_set_id NOT IN {structural_alert_set_ids});
        """
    ).format(
        assay_ids=sql.Literal(tuple(assay_ids)),
        pchembl_min_value=sql.Literal(pchembl_min_value),
        min_mw=sql.Literal(min_mw),
        max_mw=sql.Literal(max_mw),
        structural_alert_set_ids=sql.Literal(tuple(structural_alert_set_ids)),
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

    df = pl.read_csv(args.assay_tsv_file, separator="\t")
    assay_ids = list(df["assay_id"])
    mol_info = get_active_mol_info(
        cursor,
        assay_ids,
        args.pchembl_min_value,
        args.min_mw,
        args.max_mw,
        args.structural_alert_set_ids,
    )
    header = ["molregno", "chembl_id", "assay_id", "standard_type", "pchembl_value"]
    write_to_tsv(args.activities_tsv_file, mol_info, header)

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
