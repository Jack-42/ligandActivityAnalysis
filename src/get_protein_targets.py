"""
@author Jack Ringer
Date: 6/4/2025
Description:
This script takes in the name of a given protein family/class (e.g., "Kinase").
The script will return the TargetIDs (TID) of each protein belonging to the given family from ChEMBL.
"""

import argparse
import csv
import os
from typing import Optional

import psycopg2
from psycopg2 import sql
from psycopg2.extras import DictCursor


def write_to_tsv(fpath: str, rows: list, header: list[str]) -> None:
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def get_protein_class_relations(
    cursor: DictCursor,
    protein_class_id: int,
) -> list[int]:
    # Fetch protein_class_id(s) which are or fall under given ID + their parents
    query = sql.SQL(
        """
        WITH RECURSIVE protein_class_hierarchy AS (
            SELECT protein_class_id, parent_id
            FROM protein_classification
            WHERE protein_class_id = {given_id}

            UNION ALL
            
            SELECT p.protein_class_id, p.parent_id
            FROM protein_classification p
            INNER JOIN protein_class_hierarchy pch ON pch.protein_class_id = p.parent_id
        )
        SELECT * FROM protein_class_hierarchy;
        """
    ).format(given_id=sql.Literal(protein_class_id))

    cursor.execute(query)
    result = cursor.fetchall()
    return result


def get_target_info(
    cursor: DictCursor,
    protein_class_ids: list[int],
    organism: Optional[str] = None,
    target_type: Optional[str] = None,
) -> list[int]:
    query = sql.SQL(
        """
        SELECT td.tid, cc.protein_class_id
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_class cc ON tc.component_id = cc.component_id
        WHERE cc.protein_class_id IN {protein_class_ids}
        AND ({organism} IS NULL OR td.organism = {organism})
        AND ({target_type} IS NULL OR td.target_type = {target_type});
        """
    ).format(
        protein_class_ids=sql.Literal(tuple(protein_class_ids)),
        organism=sql.Literal(organism),
        target_type=sql.Literal(target_type),
    )

    cursor.execute(query)
    result = cursor.fetchall()
    return result


def parse_args():
    parser = argparse.ArgumentParser(
        description="Return the TargetIDs (TIDs) of proteins belonging to the given family, as well as protein family relation information (e.g., if the given family has subfamilies).",
        epilog="",
    )
    parser.add_argument(
        "--db_name",
        type=str,
        default=argparse.SUPPRESS,
        required=True,
        help="Name of ChEMBL PostgreSQL database",
    )
    parser.add_argument(
        "--db_host",
        type=str,
        default="localhost",
        help="Database host (default: %(default)s)",
    )
    parser.add_argument(
        "--db_user",
        type=str,
        default=argparse.SUPPRESS,
        required=True,
        help="Database user",
    )
    parser.add_argument(
        "--db_password",
        type=str,
        default=argparse.SUPPRESS,
        required=True,
        help="Database password",
    )
    parser.add_argument(
        "--db_port",
        type=int,
        default=5432,
        required=False,
        help="Database port",
    )
    parser.add_argument(
        "--protein_class_id",
        type=int,
        required=True,
        default=argparse.SUPPRESS,
        help="ID of the protein class. Should match against 'protein_class_id' column in 'protein_classification' table.",
    )
    parser.add_argument(
        "--target_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing TID and protein_class_id for each protein",
    )
    parser.add_argument(
        "--family_tsv_file",
        type=str,
        required=True,
        default=argparse.SUPPRESS,
        help="Output TSV file containing family relationship information from protein_classification table",
    )
    parser.add_argument(
        "--organism",
        type=str,
        required=False,
        default=None,
        help="(Optional) If given, will only select proteins belonging to the given organism. Make sure this arg matches the 'organism' column in the 'target_dictionary' table.",
    )
    parser.add_argument(
        "--target_type",
        type=str,
        required=False,
        default=None,
        help="(Optional) If given, will only select proteins of the given type. Make sure this arg matches the 'target_type' column in the 'target_dictionary' table.",
    )
    args = parser.parse_args()
    return args


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

    # Get class and all subclasses
    protein_class_relations = get_protein_class_relations(cursor, args.protein_class_id)
    protein_class_ids = [x[0] for x in protein_class_relations]

    # Get components of selected classes
    target_info = get_target_info(
        cursor, protein_class_ids, args.organism, args.target_type
    )

    # Write outputs
    write_to_tsv(
        args.family_tsv_file, protein_class_relations, ["protein_class_id", "parent_id"]
    )
    write_to_tsv(args.target_tsv_file, target_info, ["tid", "protein_class_id"])

    # Close connections
    cursor.close()
    connection.close()


if __name__ == "__main__":
    main()
