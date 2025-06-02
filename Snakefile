import os

from snakemake.utils import min_version

min_version("6.0")
configfile: "config.yaml"


rule all:
    input:
        config["CHEMBL_DB_FPATH"]


rule download_chembl_db:
    output:
        DB_FPATH=config["CHEMBL_DB_FPATH"]
    params:
        DB_LINK=config["CHEMBL_POSTGRES_DB_LINK"],
        OUTPUT_DIR=os.path.dirname(config["CHEMBL_DB_FPATH"])
    log:
        "logs/download_chembl_db/all.log"
    benchmark:
        "benchmark/download_chembl_db/all.tsv"
    shell:
        """
        mkdir -p {params.OUTPUT_DIR} && \
        wget -O {output.DB_FPATH} {params.DB_LINK} \
        > {log} 2>&1 
        """
