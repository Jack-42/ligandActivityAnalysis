import os
import re
from snakemake.utils import min_version

min_version("6.0")


configfile: "config.yaml"


def get_boolean_arg(key: str, param_name: str, default=None):
    return (
        f"--{param_name}" if config.get(key, default) == True else f"--no-{param_name}"
    )


CHEMBL_VERSION = re.search(r"\d+", config["CHEMBL_POSTGRES_DB_LINK"])[0]
CHEMBL_DB_TAR_FPATH = os.path.join(
    config["CHEMBL_LOCAL_DIR"], os.path.basename(config["CHEMBL_POSTGRES_DB_LINK"])
)
CHEMBL_DB_DMP_FILE = f"chembl_{CHEMBL_VERSION}_postgresql.dmp"
CHEMBL_DB_DMP_FPATH = os.path.join(
    config["CHEMBL_LOCAL_DIR"],
    f"chembl_{CHEMBL_VERSION}/chembl_{CHEMBL_VERSION}_postgresql/",
    CHEMBL_DB_DMP_FILE,
)

SETUP_DB_DONE_FILE = "logs/setup_db/done.txt"
SETUP_USER_DONE_FILE = "logs/setup_user/done.txt"


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
APT_PNG_FILE = os.path.join(FIGURES_DIR, "assays_per_target.png")
FAMILY_TREE_PNG_FILE = os.path.join(FIGURES_DIR, "protein_family_tree.png")


rule all:
    input:
        APT_PNG_FILE,
        ACTIVITIES_TSV_FILE,
        FAMILY_DETAILS_TSV_FILE,
        FAMILY_TREE_PNG_FILE,


rule download_chembl_files:
    output:
        CHEMBL_DB_TAR_FPATH=CHEMBL_DB_TAR_FPATH,
    params:
        CHEMBL_POSTGRES_DB_LINK=config["CHEMBL_POSTGRES_DB_LINK"],
        CHEMBL_LOCAL_DIR=config["CHEMBL_LOCAL_DIR"],
    log:
        "logs/download_chembl_files/all.log",
    benchmark:
        "benchmark/download_chembl_files/all.tsv"
    shell:
        """
        mkdir -p {params.CHEMBL_LOCAL_DIR} && \
        wget -O {output.CHEMBL_DB_TAR_FPATH} {params.CHEMBL_POSTGRES_DB_LINK} \
        > {log} 2>&1 
        """


rule extract_files:
    input:
        CHEMBL_DB_TAR_FPATH=CHEMBL_DB_TAR_FPATH,
    output:
        CHEMBL_DB_DMP_FPATH=CHEMBL_DB_DMP_FPATH,
    params:
        CHEMBL_LOCAL_DIR=config["CHEMBL_LOCAL_DIR"],
    log:
        "logs/extract_files/all.log",
    benchmark:
        "benchmark/extract_files/all.tsv"
    shell:
        """
        tar -xzvf {input.CHEMBL_DB_TAR_FPATH} -C {params.CHEMBL_LOCAL_DIR} \
        > {log} 2>&1 
        """


rule setup_db:
    input:
        CHEMBL_DB_DMP_FPATH=CHEMBL_DB_DMP_FPATH,
    output:
        touch(SETUP_DB_DONE_FILE),
    params:
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/setup_db/all.log",
    benchmark:
        "benchmark/setup_db/all.tsv"
    shell:
        """
        createdb -h {params.CHEMBL_DB_HOST} -p {params.CHEMBL_DB_PORT} {params.CHEMBL_DB_NAME} || echo "Database may already exist" && \
        pg_restore -O -v -h {params.CHEMBL_DB_HOST} -p {params.CHEMBL_DB_PORT} -d {params.CHEMBL_DB_NAME} {input.CHEMBL_DB_DMP_FPATH} \
        > {log} 2>&1
        """


rule setup_user:
    input:
        SETUP_DB_DONE_FILE,
    output:
        touch(SETUP_USER_DONE_FILE),
    params:
        CHEMBL_DB_HOST=config["CHEMBL_DB_HOST"],
        CHEMBL_DB_NAME=config["CHEMBL_DB_NAME"],
        CHEMBL_DB_USER=config["CHEMBL_DB_USER"],
        CHEMBL_DB_PASSWORD=config["CHEMBL_DB_PASSWORD"],
        CHEMBL_DB_PORT=config["CHEMBL_DB_PORT"],
    log:
        "logs/setup_user/all.log",
    benchmark:
        "benchmark/setup_user/all.tsv"
    shell:
        """
        psql -d {params.CHEMBL_DB_NAME} -c "CREATE ROLE {params.CHEMBL_DB_USER} WITH LOGIN PASSWORD '{params.CHEMBL_DB_PASSWORD}'" && \
        psql -d {params.CHEMBL_DB_NAME} -c "GRANT SELECT ON ALL TABLES IN SCHEMA public TO {params.CHEMBL_DB_USER}" && \
        psql -d {params.CHEMBL_DB_NAME} -c "GRANT SELECT ON ALL SEQUENCES IN SCHEMA public TO {params.CHEMBL_DB_USER}" && \
        psql -d {params.CHEMBL_DB_NAME} -c "GRANT EXECUTE ON ALL FUNCTIONS IN SCHEMA public TO {params.CHEMBL_DB_USER}" \
        > {log} 2>&1 
        """


rule get_protein_targets:
    input:
        SETUP_USER_DONE_FILE,
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
        "--interactive_mode "
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
        " > {log} 2>&1 "
