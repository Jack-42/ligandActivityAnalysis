"""
@author Jack Ringer
Date: 6/16/2025
Description:
Shared utilities for file I/O
"""

import csv
import os
import pickle


def write_to_tsv(fpath: str, rows: list, header: list[str]) -> None:
    os.makedirs(os.path.dirname(fpath), exist_ok=True)
    with open(fpath, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def save_to_pkl(obj, save_path: str) -> None:
    os.makedirs(os.path.dirname(save_path), exist_ok=True)
    with open(save_path, "wb") as f:
        pickle.dump(obj, f, protocol=5)


def load_from_pkl(save_path: str) -> object:
    with open(save_path, "rb") as f:
        obj = pickle.load(f)
    return obj
