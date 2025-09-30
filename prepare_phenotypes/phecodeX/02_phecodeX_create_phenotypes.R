library(tidyverse)
if(!require(PheWAS)) devtools::install_github("PheWAS/PheWAS")
require(PheWAS)
require(data.table)

phecode = fread("updated_phecodex_map.csv")
icd = fread("diagnosis_icd9_icd10.tsv")

data = inner_join(icd, phecode, by = c(vocabulary_id = 'vocabulary_id', concept_code = 'ICD')) |>
       group_by(person_id, phecode) |> 
       summarize(code_count = length(unique(date)))

fwrite(data, "phecode_counts.tsv", sep='\t')


ehr_inds = fread("persons_srwgs_ehr.tsv") |> filter(has_ehr_data == 1 & has_srwgs == 1)


ehr_person_ids= ehr_inds %>% select(person_id) %>% distinct() |> pull(person_id)

phe_table = createPhenotypes(data %>% transmute(person_id, vocabulary_id="phecode",phecode, code_count), 
			        translate=FALSE, 
                             min.code.count=2, 
                             add.phecode.exclusions=FALSE,
				 full.population.ids=ehr_person_ids)



fwrite(phe_table |> select(FID = person_id, IID = person_id, everything()), "phecodeX_table.tsv", sep='\t', na = "NA")

