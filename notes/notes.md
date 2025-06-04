# Notes

## ChEMBL (PostgreSQL) DB
* Visual of schema: https://ftp.ebi.ac.uk/pub/databases/
chembl/ChEMBLdb/latest/chembl_35_schema.png

### Compounds/Ligands
* Compound structures are located in the `compound_structures` table, column `molfile`. Each compound ID in column `molregno`
* Can cross-reference `compound_structures` with `molecule_dictionary` using  `molregno`
    * `molecule_dictionary` contains `chembl_id`    

### Proteins
* *Note*: The term "component" generally refers to protein subunits (see [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC3965067/))
    * Types in DB include "PROTEIN", "DNA", and "RNA"
* Relationships between "families" contained in `protein_classification` table
* Want single protein complexes: `select * from target_dictionary where target_type='SINGLE PROTEIN';`
* Get ID of kinase family (=6): `select protein_class_id from protein_classification where pref_name='Kinase';`
