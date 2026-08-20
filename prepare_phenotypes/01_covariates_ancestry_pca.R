library(tidyverse)
library(bigrquery)

PROJECT <- Sys.getenv("GOOGLE_PROJECT")
#BUCKET <- Sys.getenv("WORKSPACE_BUCKET")
CDR <- Sys.getenv("WORKSPACE_CDR")

bq_query <- function(sql) bq_table_download(bq_dataset_query(Sys.getenv("WORKSPACE_CDR"),
							     sql, 
							     billing = Sys.getenv("GOOGLE_PROJECT")),
					                     bigint = "integer64")
                        
cb_person <- bq_query("select * from cb_search_person")
cb_person$person_id = as.integer(cb_person$person_id)

ancestry <- read_tsv("/home/jupyter/workspace/vwb-aou-datasets-controlled-v9/v9/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv")

pcas <- ancestry |>
  mutate(person_id = as.integer(research_id)) |>
  select(person_id, ancestry_pred, ancestry_pred_other, pca_features) |>
  separate(pca_features, sep="[,[\\]]", into=c(NA, paste0("pc",1:16), NA))
      
covariates <- pcas |> 
  inner_join(cb_person, by = 'person_id') |>
  subset(has_ehr_data == 1) |>
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
  write_tsv("/home/jupyter/workspace/workspace-bucket/covariates.tsv")
      
covariates |>
  select(FID, IID, ancestry_pred, ancestry_pred_other) |>
  write_tsv("/home/jupyter/workspace/workspace-bucket/ancestry.tsv")
        
covariates |>
  select(FID, IID, is_male) |>
  write_tsv("/home/jupyter/workspace/workspace-bucket/is_male.tsv")
		
