"""
@author Jack Ringer
Date: 7/7/2025
Description:
Contains constant (e.g., filenames) used in multiple
locations in the project.
"""

import os


def get_ligand2cluster_fpath(cluster_dir: str, class_level: int) -> str:
    return os.path.join(cluster_dir, f"ligand2cluster-class_level={class_level}.tsv")


def get_target2cluster_fpath(cluster_dir: str, class_level: str) -> str:
    return os.path.join(cluster_dir, f"target2cluster-class_level={class_level}.tsv")


def get_cluster2sim_fpath(cluster2sim_dir: str, class_level: str) -> str:
    return os.path.join(cluster2sim_dir, f"cluster2sim_class_level={class_level}.pkl")
