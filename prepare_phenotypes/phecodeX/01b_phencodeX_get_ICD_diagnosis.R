library(tidyverse)
library(bigrquery)

GOOGLE_PROJECT <- Sys.getenv("GOOGLE_PROJECT")
WORKSPACE_BUCKET <- "gs://workspace-bucket-wb-radiant-cabbage-3726/users/florian/projects/phecodex_phenotypes"
WORKSPACE_CDR <- Sys.getenv("WORKSPACE_CDR")

bq_export_tsv <- function(sql, gcs_prefix, local_out) {
  # 1. Run query into a temporary BigQuery table
  tb <- bq_dataset_query(
    WORKSPACE_CDR,
    sql,
    billing = GOOGLE_PROJECT
  )
  
  # 2. Export directly from BigQuery to your All of Us bucket (sharded if > 1 GB)
  gcs_uri <- glue::glue("{WORKSPACE_BUCKET}/bq_exports/{gcs_prefix}-*.tsv")
  bq_table_save(
    x = tb,
    destination_uris = gcs_uri,
    destination_format = "CSV",
    field_delimiter = "\t",
    print_header = TRUE
  )
  
  # 3. Combine shards cleanly into a single local TSV (keeping header only once)
  tmp_dir <- tempfile()
  dir.create(tmp_dir)
  system(glue::glue("gsutil -m cp '{gcs_uri}' {tmp_dir}/"))
  
  files <- sort(list.files(tmp_dir, full.names = TRUE))
  file.copy(files[1], local_out, overwrite = TRUE)
  if (length(files) > 1) {
    for (f in files[-1]) {
      system(glue::glue("tail -n +2 {shQuote(f)} >> {shQuote(local_out)}"))
    }
  }
  unlink(tmp_dir, recursive = TRUE)
}

sql1 <- str_glue("
  select distinct person_id, vocabulary_id, concept_code, condition_start_date as date
  from `{WORKSPACE_CDR}.condition_occurrence`
  join `{WORKSPACE_CDR}.concept` on (condition_source_concept_id = concept_id)
  where vocabulary_id in ('ICD9CM','ICD10CM')
  union distinct
  select distinct person_id, vocabulary_id, concept_code, observation_date as date
  from `{WORKSPACE_CDR}.observation`
  join `{WORKSPACE_CDR}.concept` on (observation_source_concept_id = concept_id)
  where vocabulary_id in ('ICD9CM','ICD10CM')
  union distinct
  select distinct person_id, vocabulary_id, concept_code, procedure_date as date
  from `{WORKSPACE_CDR}.procedure_occurrence`
  join `{WORKSPACE_CDR}.concept` on (procedure_source_concept_id = concept_id)
  where vocabulary_id in ('ICD9CM','ICD10CM')
  union distinct
  select distinct person_id, vocabulary_id, concept_code, measurement_date as date
  from `{WORKSPACE_CDR}.measurement`
  join `{WORKSPACE_CDR}.concept` on (measurement_source_concept_id = concept_id)
  where vocabulary_id in ('ICD9CM','ICD10CM')
")

bq_export_tsv(sql1, "diagnosis_icd9_icd10", "diagnosis_icd9_icd10.csv")

sql2 <- str_glue("
  select person_id, has_whole_genome_variant as has_srwgs, has_ehr_data
  from `{WORKSPACE_CDR}.cb_search_person`
")

bq_export_tsv(sql2, "persons_srwgs_ehr", "persons_srwgs_ehr.csv")

