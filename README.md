# ligandActivityAnalysis

Project for Jack Ringer as part of the CFDE 2025 Internship @UNM

## PostgreSQL Setup

**Note:** This workflow has been tested with PostGreSQL version `psql (PostgreSQL) 16.9 (Ubuntu 16.9-0ubuntu0.24.04.1)`

If you do not already have PostGreSQL installed, you can follow the instructions [here](https://www.postgresql.org/download/).

After installing PostGreSQL, you need to make your user a superuser prior to DB setup:

1. Switch to postgres user: `(base) <username>@<computer>:~$ sudo -u postgres psql`

2) Make yourself a superuser: `CREATE ROLE "<username>" WITH SUPERUSER PASSWORD '<password>'`
3) Enable login: `ALTER ROLE "<username>" WITH LOGIN;`

## Running Snakemake Workflow

1. Make sure to setup PostgreSQL as above
2. Install conda environment: `conda env create -f environment.yaml`
3. Modify `config.yaml` to your desired params
4. Run snakemake: `snakemake --cores 1`
   - If prompted for `Password:`, enter the `<password>` of your SUPERUSER
