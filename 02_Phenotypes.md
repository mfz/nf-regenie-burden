## Phenotypes

### PhecodeX phenotypes

Run the scripts in `prepare_phenotypes/phecodeX` to generate PhecodeX phenotypes.
This creates a file phecodeX_table.tsv with columns FID, IID,  and one column per phecode.



### OMOP case-control phenotypes

Run the scripts  

- `prepare_phenotypes/OMOP_cc/01_OMOP_get_data.R`
- `prepare_phenotypes/OMOP_cc/02_OMOP_create_database.R`

to create a database OMOP.db with OMOP concepts.

This database can then be used to create case-control phenotypes based on condition and procedure concepts.

See `prepare_phenotypes/OMOP_cc/03_OMOPdb_query_API.ipynb` for an example of how to create case-control phenotypes.







