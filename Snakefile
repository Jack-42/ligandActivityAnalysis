import os
import re
from snakemake.utils import min_version

min_version("6.0")


configfile: "config.yaml"


def get_boolean_arg(key: str, param_name: str, default=None):
    return (
        f"--{param_name}" if config.get(key, default) == True else f"--no-{param_name}"
    )


# want to make it easy to run workflow with different organism, target_type, etc
# all other data depends on selected targets
SUBDIR_NAME = "".join(
    f'{config["ORGANISM"]}-{config["TARGET_TYPE"]}-{config["PROTEIN_CLASS_ID"]}'.split()
)
DATA_SUBDIR = os.path.join(config["DATA_DIR"], SUBDIR_NAME)
FIGURES_DIR = os.path.join(DATA_SUBDIR, "figures")
TARGET_DATA_DIR = os.path.join(DATA_SUBDIR, "targets")
TARGET_TSV_FILE = os.path.join(TARGET_DATA_DIR, "target_info.tsv")
FAMILY_TSV_FILE = os.path.join(TARGET_DATA_DIR, "family_info.tsv")
FAMILY_DETAILS_TSV_FILE = os.path.join(TARGET_DATA_DIR, "family_addl_info.tsv")
ASSAYS_TSV_FILE = os.path.join(DATA_SUBDIR, "assays", "assay_info.tsv")
ACTIVITIES_TSV_FILE = os.path.join(DATA_SUBDIR, "activities", "activity_info.tsv")
STRUCTURE_SUBDIR = os.path.join(DATA_SUBDIR, "compound_structures")
COMPOUND_STRUCTURES_TSV_FILE = os.path.join(STRUCTURE_SUBDIR, "structure_info.tsv")
COMPOUND_FINGERPRINTS_PKL_FILE = os.path.join(STRUCTURE_SUBDIR, "fingerprints.pkl")
APT_PNG_FILE = os.path.join(FIGURES_DIR, "assays_per_target.png")
FAMILY_TREE_PNG_FILE = os.path.join(FIGURES_DIR, "protein_family_tree.png")

# ligands clustered by classification of their target(s)
# as well as targets "clustered" by family classification
CLUSTER_DIR = os.path.join(DATA_SUBDIR, "family_clusters")
LIGAND2TID_TSV_FILE = os.path.join(CLUSTER_DIR, "ligand2tid.tsv")

SIMILARITY_NPY_FILE = os.path.join(
    STRUCTURE_SUBDIR, "similarity_lower_triangular_matrix.npy"
)
SIMILARITY_ID_PKL_FILE = os.path.join(STRUCTURE_SUBDIR, "similarity_id_list.pkl")
CLUSTER2SIM_DIR = os.path.join(DATA_SUBDIR, "processed_similarity_cluster_data")


rule all:
    input:
        APT_PNG_FILE,
        FAMILY_DETAILS_TSV_FILE,
        FAMILY_TREE_PNG_FILE,
        CLUSTER2SIM_DIR,


rule get_protein_targets:
    output:
        target_tsv_file=TARGET_TSV_FILE,
        family_tsv_file=FAMILY_TSV_FILE,
    params:
        ORGANISM=config["ORGANISM"],
        TARGET_TYPE=config["TARGET_TYPE"],
        PROTEIN_CLASS_ID=config["PROTEIN_CLASS_ID"],
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/get_protein_targets/all.log",
    benchmark:
        "benchmark/get_protein_targets/all.tsv"
    shell:
        "python src/get_protein_targets.py "
        "--db_name '{params.CHEMBL_DB_NAME}' "
        "--db_host '{params.CHEMBL_DB_HOST}' "
        "--db_user '{params.CHEMBL_DB_USER}' "
        "--db_password '{params.CHEMBL_DB_PASSWORD}' "
        "--db_port {params.CHEMBL_DB_PORT} "
        "--protein_class_id {params.PROTEIN_CLASS_ID} "
        "--target_type '{params.TARGET_TYPE}' "
        "--organism '{params.ORGANISM}' "
        "--target_tsv_file '{output.target_tsv_file}' "
        "--family_tsv_file '{output.family_tsv_file}' "
        " > {log} 2>&1 "


rule get_addl_family_info:
    input:
        family_tsv_file=FAMILY_TSV_FILE,
    output:
        family_details_tsv_file=FAMILY_DETAILS_TSV_FILE,
    params:
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/get_addl_family_info/all.log",
    benchmark:
        "benchmark/get_addl_family_info/all.tsv"
    shell:
        "python src/get_addl_family_info.py "
        "--db_name '{params.CHEMBL_DB_NAME}' "
        "--db_host '{params.CHEMBL_DB_HOST}' "
        "--db_user '{params.CHEMBL_DB_USER}' "
        "--db_password '{params.CHEMBL_DB_PASSWORD}' "
        "--db_port {params.CHEMBL_DB_PORT} "
        "--family_tsv_file '{input.family_tsv_file}' "
        "--family_details_tsv_file '{output.family_details_tsv_file}' "
        " > {log} 2>&1 "


rule get_relevant_assays:
    input:
        target_tsv_file=TARGET_TSV_FILE,
    output:
        assay_tsv_file=ASSAYS_TSV_FILE,
    params:
        ASSAY_TYPE=config["ASSAY_TYPE"],
        DOC_TYPE=config["DOC_TYPE"],
        CONFIDENCE_SCORE=config["CONFIDENCE_SCORE"],
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
        EXCLUDE_VARIANTS=get_boolean_arg("EXCLUDE_VARIANTS", "exclude_variants", False),
    log:
        "logs/get_relevant_assays/all.log",
    benchmark:
        "benchmark/get_relevant_assays/all.tsv"
    shell:
        "python src/get_relevant_assays.py "
        "--db_name '{params.CHEMBL_DB_NAME}' "
        "--db_host '{params.CHEMBL_DB_HOST}' "
        "--db_user '{params.CHEMBL_DB_USER}' "
        "--db_password '{params.CHEMBL_DB_PASSWORD}' "
        "--db_port {params.CHEMBL_DB_PORT} "
        "--target_tsv_file '{input.target_tsv_file}' "
        "--assay_tsv_file '{output.assay_tsv_file}' "
        "--assay_type '{params.ASSAY_TYPE}' "
        "--confidence_score {params.CONFIDENCE_SCORE} "
        "--doc_type '{params.DOC_TYPE}' "
        "{params.EXCLUDE_VARIANTS} "
        " > {log} 2>&1 "


rule plot_assays_per_target:
    input:
        assay_tsv_file=ASSAYS_TSV_FILE,
    output:
        out_path=APT_PNG_FILE,
    log:
        "logs/plot_assays_per_target/all.log",
    benchmark:
        "benchmark/plot_assays_per_target/all.tsv"
    shell:
        "python src/plotting/plot_assays_per_target.py "
        "--assay_tsv_file '{input.assay_tsv_file}' "
        "--out_path '{output.out_path}' "
        " > {log} 2>&1 "


rule visualize_protein_family_tree:
    input:
        family_details_tsv_file=FAMILY_DETAILS_TSV_FILE,
    output:
        out_path=FAMILY_TREE_PNG_FILE,
    log:
        "logs/visualize_protein_family_tree/all.log",
    benchmark:
        "benchmark/visualize_protein_family_tree/all.tsv"
    shell:
        "python src/visualize_protein_family_tree.py "
        "--family_details_tsv_file '{input.family_details_tsv_file}' "
        "--out_path '{output.out_path}' "
        " > {log} 2>&1 "


rule get_active_compounds:
    input:
        assay_tsv_file=ASSAYS_TSV_FILE,
    output:
        activities_tsv_file=ACTIVITIES_TSV_FILE,
    params:
        PCHEMBL_MIN_VALUE=config["PCHEMBL_MIN_VALUE"],
        MIN_MW=config["MIN_MW"],
        MAX_MW=config["MAX_MW"],
        STRUCTURAL_ALERT_SET_IDS=config["STRUCTURAL_ALERT_SET_IDS"],
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/get_active_compounds/all.log",
    benchmark:
        "benchmark/get_active_compounds/all.tsv"
    shell:
        "python src/get_active_compounds.py "
        "--db_name '{params.CHEMBL_DB_NAME}' "
        "--db_host '{params.CHEMBL_DB_HOST}' "
        "--db_user '{params.CHEMBL_DB_USER}' "
        "--db_password '{params.CHEMBL_DB_PASSWORD}' "
        "--db_port {params.CHEMBL_DB_PORT} "
        "--assay_tsv_file '{input.assay_tsv_file}' "
        "--pchembl_min_value {params.PCHEMBL_MIN_VALUE} "
        "--activities_tsv_file '{output.activities_tsv_file}' "
        "--min_mw {params.MIN_MW} "
        "--max_mw {params.MAX_MW} "
        "--structural_alert_set_ids {params.STRUCTURAL_ALERT_SET_IDS} "
        " > {log} 2>&1 "


rule get_compound_structures:
    input:
        activities_tsv_file=ACTIVITIES_TSV_FILE,
    output:
        compound_structures_tsv_file=COMPOUND_STRUCTURES_TSV_FILE,
    params:
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/get_compound_structures/all.log",
    benchmark:
        "benchmark/get_compound_structures/all.tsv"
    shell:
        "python src/get_compound_structures.py "
        "--db_name '{params.CHEMBL_DB_NAME}' "
        "--db_host '{params.CHEMBL_DB_HOST}' "
        "--db_user '{params.CHEMBL_DB_USER}' "
        "--db_password '{params.CHEMBL_DB_PASSWORD}' "
        "--db_port {params.CHEMBL_DB_PORT} "
        "--activities_tsv_file '{input.activities_tsv_file}' "
        "--compound_structures_tsv_file '{output.compound_structures_tsv_file}' "
        " > {log} 2>&1 "


rule generate_fingerprints:
    input:
        compound_structures_tsv_file=COMPOUND_STRUCTURES_TSV_FILE,
    output:
        fingerprints_pkl_file=COMPOUND_FINGERPRINTS_PKL_FILE,
    log:
        "logs/generate_fingerprints/all.log",
    benchmark:
        "benchmark/generate_fingerprints/all.tsv"
    shell:
        "python src/generate_fingerprints.py "
        "--compound_structures_tsv_file '{input.compound_structures_tsv_file}' "
        "--fingerprints_pkl_file '{output.fingerprints_pkl_file}' "
        " > {log} 2>&1 "


rule assign_family_clusters:
    input:
        activities_tsv_file=ACTIVITIES_TSV_FILE,
        assay_tsv_file=ASSAYS_TSV_FILE,
        target_tsv_file=TARGET_TSV_FILE,
        family_details_tsv_file=FAMILY_DETAILS_TSV_FILE,
    output:
        cluster_dir=directory(CLUSTER_DIR),
        ligand2tid_tsv_file=LIGAND2TID_TSV_FILE,
    params:
        min_class_level=config["MIN_CLASS_LEVEL"],
        max_class_level=config["MAX_CLASS_LEVEL"],
    log:
        "logs/assign_family_clusters/all.log",
    benchmark:
        "benchmark/assign_family_clusters/all.tsv"
    shell:
        "python src/assign_family_clusters.py "
        "--activities_tsv_file '{input.activities_tsv_file}' "
        "--assay_tsv_file '{input.assay_tsv_file}' "
        "--target_tsv_file '{input.target_tsv_file}' "
        "--family_details_tsv_file '{input.family_details_tsv_file}' "
        "--cluster_dir '{output.cluster_dir}' "
        "--ligand2tid_tsv_file '{output.ligand2tid_tsv_file}' "
        "--min_class_level {params.min_class_level} "
        "--max_class_level {params.max_class_level} "
        " > {log} 2>&1 "


rule calculate_fp_similarity:
    input:
        fingerprints_pkl_file=COMPOUND_FINGERPRINTS_PKL_FILE,
    output:
        similarity_npy_file=SIMILARITY_NPY_FILE,
        similarity_id_pkl_file=SIMILARITY_ID_PKL_FILE,
    params:
        max_n_compounds=config["MAX_N_COMPOUNDS"],
    log:
        "logs/calculate_fp_similarity/all.log",
    benchmark:
        "benchmark/calculate_fp_similarity/all.tsv"
    shell:
        "python src/calculate_fp_similarity.py "
        "--fingerprints_pkl_file '{input.fingerprints_pkl_file}' "
        "--similarity_npy_file '{output.similarity_npy_file}' "
        "--similarity_id_pkl_file '{output.similarity_id_pkl_file}' "
        "--max_n_compounds {params.max_n_compounds} "
        " > {log} 2>&1 "


# note: can include fingerprints file as input to validate script logic (at the cost of extra compute)
# (I've already validated this with a large matrix, so it should be ok)
rule process_similarity_by_cluster:
    input:
        similarity_npy_file=SIMILARITY_NPY_FILE,
        similarity_id_pkl_file=SIMILARITY_ID_PKL_FILE,
        cluster_dir=CLUSTER_DIR,
    output:
        cluster2sim_dir=directory(CLUSTER2SIM_DIR),
    params:
        min_class_level=config["MIN_CLASS_LEVEL"],
        max_class_level=config["MAX_CLASS_LEVEL"],
    log:
        "logs/process_similarity_by_cluster/all.log",
    benchmark:
        "benchmark/process_similarity_by_cluster/all.tsv"
    shell:
        "python src/process_similarity_by_cluster.py "
        "--similarity_npy_file '{input.similarity_npy_file}' "
        "--similarity_id_pkl_file '{input.similarity_id_pkl_file}' "
        "--cluster_dir '{input.cluster_dir}' "
        "--cluster2sim_dir '{output.cluster2sim_dir}' "
        "--min_class_level {params.min_class_level} "
        "--max_class_level {params.max_class_level} "
        " > {log} 2>&1 "
