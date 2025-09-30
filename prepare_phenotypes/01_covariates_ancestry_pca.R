library(tidyverse)
library(bigrquery)

PROJECT <- Sys.getenv("GOOGLE_PROJECT")
BUCKET <- Sys.getenv("WORKSPACE_BUCKET")
CDR <- Sys.getenv("WORKSPACE_CDR")

bq_query <- function(sql) bq_table_download(bq_dataset_query(Sys.getenv("WORKSPACE_CDR"),
							     sql, 
							     billing = Sys.getenv("GOOGLE_PROJECT")),
					                     bigint = "integer64")
                        
cb_person <- bq_query("select * from cb_search_person")

system(str_glue("gsutil -u {PROJECT} cp gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv ."), intern = T)

ancestry <- read_tsv("ancestry_preds.tsv")

pcas <- ancestry |>
  mutate(person_id = as.integer(research_id)) |>
  select(person_id, ancestry_pred, ancestry_pred_other, pca_features) |>
  separate(pca_features, sep="[,[\\]]", into=c(NA, paste0("pc",1:16), NA))
      
covariates <- pcas |> 
  inner_join(cb_person, by = 'person_id') |>
  subset(has_whole_genome_variant == 1 & has_array_data == 1) |>
  select(-starts_with("has_fitbit")) |>
  subset(sex_at_birth %in% c('Male', 'Female')) |>
  mutate(sex = sex_at_birth,
	 is_male = ifelse(sex_at_birth == 'Male', 1, 0),
	 age = age_at_cdr,
         age2 = age^2) |>
  mutate(across(pc1:pc16, as.numeric)) |>
  select(FID = person_id, IID = person_id, age, age2, sex, is_male, pc1:pc16, everything())
    
covariates |>
  select(FID, IID, is_male, age, age2, pc1:pc16) |>
  write_tsv("covariates.tsv")
      
covariates |>
  select(FID, IID, ancestry_pred, ancestry_pred_other) |>
  write_tsv("ancestry.tsv")
        
covariates |>
  select(FID, IID, is_male) |>
  write_tsv("is_male.tsv")
		
system(str_glue("gsutil cp covariates.tsv {BUCKET}/covariates.tsv"))
system(str_glue("gsutil cp ancestry.tsv {BUCKET}/ancestry.tsv"))
system(str_glue("gsutil cp is_male.tsv {BUCKET}/is_male.tsv"))
