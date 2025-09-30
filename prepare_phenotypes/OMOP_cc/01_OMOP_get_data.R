library(tidyverse)
library(bigrquery)
library(data.table)

GOOGLE_PROJECT <- Sys.getenv("GOOGLE_PROJECT")
WORKSPACE_BUCKET <- Sys.getenv("WORKSPACE_BUCKET")
WORKSPACE_CDR <- Sys.getenv("WORKSPACE_CDR")


bq_query <- function(sql) bq_table_download(bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), 
							     sql, 									       						billing = Sys.getenv("GOOGLE_PROJECT")),
						  	     bigint = "integer64")

# condition_occurrence
sql = "SELECT
    co.person_id,
    co.condition_concept_id,
    co.condition_start_date
    FROM condition_occurrence co
    JOIN concept c
        ON co.condition_concept_id = c.concept_id
    WHERE c.domain_id = 'Condition'
      AND c.standard_concept = 'S'"

df = bq_query(sql)

if (!dir.exists("OMOP") dir.create("OMOP")
fwrite(df, "OMOP/condition_occurrence.tsv", sep = "\t")

# condition_descendants
descendants_sql = "SELECT
    ca.ancestor_concept_id,
    ca.descendant_concept_id,
    d.concept_name AS descendant_name
    FROM concept_ancestor ca
    JOIN concept a 
        ON ca.ancestor_concept_id = a.concept_id
    JOIN concept d 
        ON ca.descendant_concept_id = d.concept_id
    WHERE a.domain_id = 'Condition'
      AND a.standard_concept = 'S'
      AND d.standard_concept = 'S'
      AND d.domain_id = 'Condition'"

descendants = bq_query(descendants_sql)

fwrite(descendants, "OMOP/condition_descendants.tsv", sep = "\t")

# icd_to_omop
icd_to_omop_sql = "SELECT
      icd.concept_id AS icd_concept_id,
      icd.concept_code AS icd_code,
      icd.concept_name AS icd_name,
      icd.vocabulary_id AS icd_vocabulary,
      omop.concept_id AS omop_concept_id,
      omop.vocabulary_id AS source_vocabulary,  
      omop.concept_code AS source_code,
      omop.concept_name AS source_name
      FROM concept icd
      JOIN concept_relationship cr
	ON icd.concept_id = cr.concept_id_1
	AND cr.relationship_id = 'Maps to'
      JOIN concept omop
      ON cr.concept_id_2 = omop.concept_id
      WHERE icd.vocabulary_id IN ('ICD9CM', 'ICD10CM')
	AND omop.standard_concept = 'S'
	AND omop.domain_id = 'Condition'"

icd_to_omop = bq_query(icd_to_omop_sql)

fwrite(icd_to_omop, "OMOP/icd_to_omop.tsv", sep = "\t")

# demographics
sql = "select person_id, 
    has_whole_genome_variant as has_srwgs, 
    has_ehr_data, 
    dob, 
    age_at_consent,
    sex_at_birth 
    from cb_search_person"

demographics <- bq_query(sql)

fwrite(demographics, "OMOP/demographics.tsv", sep = "\t")

# ancestry
system(str_glue("gsutil -u {GOOGLE_PROJECT} cp gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv OMOP/ancestry.tsv"), intern = T)


# procedure_occurrence
sql = "SELECT
    po.person_id,
    po.procedure_concept_id,
    po.procedure_date
    FROM procedure_occurrence po
    JOIN concept c
        ON po.procedure_concept_id = c.concept_id
    WHERE c.domain_id = 'Procedure'
      AND c.standard_concept = 'S'"

po = bq_query(sql)

fwrite(po, "OMOP/procedure_occurrence.tsv", sep='\t')

# procedure_descendants
descendants_sql = "SELECT
    ca.ancestor_concept_id,
    ca.descendant_concept_id,
    d.concept_name AS descendant_name
    FROM concept_ancestor ca
    JOIN concept a 
        ON ca.ancestor_concept_id = a.concept_id
    JOIN concept d 
        ON ca.descendant_concept_id = d.concept_id
    WHERE a.domain_id = 'Procedure'
      AND a.standard_concept = 'S'
      AND d.standard_concept = 'S'
      AND d.domain_id = 'Procedure'"

descendants = bq_query(descendants_sql)

fwrite(descendants, "OMOP/procedure_descendants.tsv", sep = '\t')

