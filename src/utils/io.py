"""
@author Jack Ringer
Date: 6/16/2025
Description:
Shared utilities for file I/O
"""

import csv
import os


def write_to_tsv(fpath: str, rows: list, header: list[str]) -> None:
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
