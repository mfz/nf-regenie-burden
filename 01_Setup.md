## Installation
### Setting up nextflow

We need a recent version of nextflow which we install under `~` (HOME directory).

```sh
curl -s https://get.sdkman.io | bash
sdk install java 17.0.10-tem
curl -s https://get.nextflow.io | bash
```


### Installing the nextflow pipeline

```sh
git clone git@github.com:mfz/nf-regenie-burden.git
```

### Creating split ACAF bgen files

**NOTE**: This takes about 4 TB of space in the bucket, which costs around 80$/month. It is, therefore, recommended to do this only in one workspace bucket. So check if this already exists!

Google Cloud Batch does not have access to the AllOfUS genotype buckets. So we need to copy them to our own workspace buckets. We also want to split them into chunks of 200 000 variants, such that the jobs can be run within an hour or so. This allows us to use the cheaper spot instances.

Here we split the ACAF bgen files. ACAF stands for Allele Count Allele Frequency. These bgen files only include variants that have Allele Count > 100 or Allele Freqeuncy > 1%.

```sh
# confirm that this works!
# we might need to copy the original bgen files to our workspace bucket first
# to be able to access them! 
cd nf-regenie-burden/utils
~/nextflow run splitbgen.nf -profile gcb -with-report splitbgen.html
```

This installs the split ACAF bgen files into ${WORKSPACE_BUCKET}/split_bgen.
One can check this using `gsutil ls ${WORKSPACE_BUCKET}/split_bgen`.

### Create exome-specific bgen files

AllOfUS also provides exome bgen files, which are not filtered based on allele counts/frequency. There is some overap with the ACAF files, though. Here we create exome-specific bgen files that only contain those variants in the exome regions that are not already included in the ACAF files. We also split them into chunks of 200 000 variants.

```sh
gsutil cp  gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/bgen/*.bgi  ${WORKSPACE_BUCKET}/acaf_bgi

gsutil cp  gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/bgen/* ${WORKSPACE_BUCKET}/exome_bgen

cd nf-regenie-burden/utils
~/nextflow run exome_specific.nf -profile gcb -with-report exome_specific.html

gsutil cp ${WORKSPACE_BUCKET}/exome_only/*.bgen ${WORKSPACE_BUCKET}/split_bgen

```

**NOTE**: If you want to run GWAS on exomes only, you should copy the original exome bgen files to your workspace bucket and use those.

```sh
# as already done above
gsutil cp  gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/exome/bgen/* ${WORKSPACE_BUCKET}/exome_bgen
```


### Get covariates, ancestry and PCAs

_phenotype_setup/01_covariates_ancestry_pca.R_

```R
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
```


### Split PLINK files by ancestry and sex

Adjust params in `splitancestry.nf` before running

```sh
# we probably need to copy the PLINK data to our own workspace
# in order to access it ...
cd nf-regenie-burden/utils
~/nextflow run splitancestry.nf -profile gcb -with-report splitancestry.html 
```


### Create burden annotation files

see BID-2321








