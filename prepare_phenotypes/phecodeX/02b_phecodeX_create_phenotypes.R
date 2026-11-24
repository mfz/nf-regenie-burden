require(tidyverse)
require(data.table)

phecode = fread("updated_phecodex_map.csv")
icd = fread("diagnosis_icd9_icd10.csv")

genotyped = fread("/home/jupyter/workspace/workspace-bucket/data/exome_only/chr1_exome_only.sample") |> select(person_id = ID_1) # 535k

ehr_inds = fread("persons_srwgs_ehr.csv") |> filter(has_ehr_data == 1) |> select(person_id) # 481k
ancestry = fread("/home/jupyter/workspace/workspace-bucket/data/ancestry.tsv") |> select(person_id=IID, ancestry=ancestry_pred) # 535k only !!!!
covariates = fread("/home/jupyter/workspace/workspace-bucket/data/covariates.tsv") |> select(person_id=IID, sex_male = is_male)

ancestry |> group_by(ancestry) |> count()


universe <- ehr_inds |> 
    inner_join(ancestry) |>
    inner_join(covariates) |>
    select(person_id) |> distinct()


data = inner_join(icd, phecode, by = c(vocabulary_id = 'vocabulary_id', concept_code = 'ICD')) |>
       distinct(person_id, phecode) |>
       inner_join(universe)

counts = data |>
       inner_join(ancestry) |>
       inner_join(covariates) |>
       group_by(ancestry, phecode) |>
       summarise(ncases = n(), .groups="drop_last")

counts.sex = data |>
       inner_join(ancestry) |>
       inner_join(covariates) |>
       group_by(ancestry, sex_male, phecode) |>
       summarise(ncases = n(), .groups="drop_last")

fwrite(counts, file="phecodex_counts_by_ancestry.tsv", sep='\t')
fwrite(counts.sex, file="phecodex_counts_by_ancestry_sex.tsv", sep='\t')

data <- data |> 
       mutate(value = 1) |>
       pivot_wider(names_from = phecode, values_from = value, values_fill = 0)

# Ensure all subjects in universe included
data <- universe %>%
  left_join(data, by = "person_id") %>%
  mutate(across(-person_id, ~ replace_na(.x, 0)))

fwrite(data |> select(FID=person_id, IID=person_id, everything()), file="phecodex_matrix.tsv", sep='\t')



