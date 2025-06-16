"""
@author Jack Ringer
Date: 6/16/2025
Description:
Shared arguments.
"""

import argparse


def add_chembl_db_args(parser: argparse.ArgumentParser) -> None:
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
