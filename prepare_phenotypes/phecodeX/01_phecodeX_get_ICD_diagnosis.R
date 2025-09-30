library(tidyverse)
library(bigrquery)

GOOGLE_PROJECT <- Sys.getenv("GOOGLE_PROJECT")
WORKSPACE_BUCKET <- Sys.getenv("WORKSPACE_BUCKET")
WORKSPACE_CDR <- Sys.getenv("WORKSPACE_CDR")


bq_query <- function(sql) bq_table_download(bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), 
			                   sql, 
	                                   billing = Sys.getenv("GOOGLE_PROJECT")),
		                           bigint = "integer64")

sql = str_glue("select * from ( \
		  select distinct person_id, vocabulary_id, concept_code, condition_start_date as date \
		  from {WORKSPACE_CDR}.condition_occurrence \
		  join {WORKSPACE_CDR}.concept on (condition_source_concept_id=concept_id) \
		  where vocabulary_id in ('ICD9CM','ICD10CM') \
		Union distinct \
		  select distinct person_id, vocabulary_id, concept_code, observation_date as date \
		  from {WORKSPACE_CDR}.observation \
		  join {WORKSPACE_CDR}.concept on (observation_source_concept_id=concept_id) \
		  where vocabulary_id in ('ICD9CM','ICD10CM') \
		union distinct \
		  select distinct person_id, vocabulary_id, concept_code, procedure_date as date \
		  from {WORKSPACE_CDR}.procedure_occurrence \
		  join {WORKSPACE_CDR}.concept on (procedure_source_concept_id=concept_id) \
		  where vocabulary_id in ('ICD9CM','ICD10CM') \
		union distinct \
		  select distinct person_id, vocabulary_id, concept_code, measurement_date as date \
		  from {WORKSPACE_CDR}.measurement \
		  join {WORKSPACE_CDR}.concept on (measurement_source_concept_id=concept_id) \
		  where vocabulary_id in ('ICD9CM','ICD10CM'))")

data = bq_query(sql)

require(data.table)
fwrite(data, "diagnosis_icd9_icd10.tsv", sep = "\t", quote = FALSE)

sql = "select person_id, has_whole_genome_variant as has_srwgs, has_ehr_data from cb_search_person"

persons = bq_query(sql)
fwrite(persons, "persons_srwgs_ehr.tsv", sep = "\t", quote = FALSE)
	        
