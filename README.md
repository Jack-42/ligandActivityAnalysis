# ligandActivityAnalysis

Project for Jack Ringer as part of the CFDE 2025 Internship @UNM.

Although this project has been used to compare the computational similarities of ligands active within different protein kinase groups, the Snakemake workflow could be adapted to explore other protein families or to assist with other projects using ChEMBL data.

## About (Abstract)
This work investigates whether there is a relationship between the 2D computational similarity of ligands and their activity within specific human protein kinase groups. Using data from ChEMBL binding assays, the distribution of pairwise Tanimoto similarity coefficients values computed between all protein kinase ligands was compared against the distribution of pairwise similarity values computed with respect to ligands active within a specific protein kinase group. Except for the CK1 group, no significant group-specific differences were found. These results suggest there is limited utility of 2D similarity metrics for identifying ligand selectivity across a majority of protein kinase groups. However, given the many confounders that exist when performing large scale computational analyses of ChEMBL bioassay data these results are not definitive, and limitations as well as potential follow-ups are discussed in detail. 

![alt text](https://github.com/Jack-42/ligandActivityAnalysis/blob/main/docs/figures/violin_plot.png)


## Running Snakemake Workflow
1. Setup PostgreSQL on your system by following [these instructions](#postgresql-setup)
2. Install Conda environment: `conda env create -f environment.yaml`
3. Modify `config.yaml` to your desired params
4. (Optional): Run `export PGPASSWORD=<your_password>` - avoids password prompts during the DB build process
5. If you don't have it installed already, set up the ChEMBL PostgreSQL DB: `snakemake --snakefile Snakefile_DB_Setup --cores 1`
6. Run main Snakemake workflow: `snakemake --cores 1`
   - If prompted for `Password:`, enter the `<password>` of your SUPERUSER

### Workflow diagram
Output of `snakemake --forceall --rulegraph | dot -Tpng > rulegraph.png`:

![alt text](https://github.com/Jack-42/ligandActivityAnalysis/blob/main/docs/figures/rulegraph.png)


## PostgreSQL Setup

**Note:** This workflow has been tested with PostgreSQL version `psql (PostgreSQL) 16.9 (Ubuntu 16.9-0ubuntu0.24.04.1)`

If you do not already have PostgreSQL installed, you can follow the instructions [here](https://www.postgresql.org/download/).

After installing PostgreSQL, you need to make your user a superuser prior to DB setup:

1. Switch to postgres user: `(base) <username>@<computer>:~$ sudo -u postgres psql`

2) Make yourself a superuser: `CREATE ROLE "<username>" WITH SUPERUSER PASSWORD '<password>'`
3) Enable login: `ALTER ROLE "<username>" WITH LOGIN;`
