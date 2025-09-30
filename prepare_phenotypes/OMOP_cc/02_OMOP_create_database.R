require(RSQLite)
require(data.table)
require(dplyr)

con <- dbConnect(SQLite(), "OMOP/OMOP.db")

co <- fread("OMOP/condition_occurrence.tsv")
dbWriteTable(con, "condition_occurrence", co, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_co_cci ON condition_occurrence(condition_concept_id)")
co <- NULL

po <- fread("OMOP/procedure_occurrence.tsv")
dbWriteTable(con, "procedure_occurrence", po, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_po_pci ON procedure_occurrence(procedure_concept_id)")
po <- NULL

cd <- fread("OMOP/condition_descendants.tsv")
dbWriteTable(con, "condition_descendants", cd, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_cd_aci ON condition_descendants(ancestor_concept_id)")
cd <- NULL

pd <- fread("OMOP/procedure_descendants.tsv")
dbWriteTable(con, "procedure_descendants", pd, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_pd_aci ON procedure_descendants(ancestor_concept_id)")
pd <- NULL

icd2omop <- fread("OMOP/icd_to_omop.tsv")
dbWriteTable(con, "icd2omop", icd2omop, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_i2o_icdcode_voc ON icd2omop(icd_code, icd_vocabulary)")
icd2omop <- NULL

demo <- fread("OMOP/demographics.tsv")
ancestry <- fread("OMOP/ancestry.tsv") |> select(person_id = research_id, ancestry_pred)
demo <- left_join(demo, ancestry)
dbWriteTable(con, "demographics", demo, overwrite = TRUE)
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_demographics_pi ON demographics(person_id)")
demo <- NULL

dbDisconnect(con)


