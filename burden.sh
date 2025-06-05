
genotype_array="../arrays/arrays_EUR_qcfiltered"
phenotype_file="../phecodeX_test.tsv"
covariates_file="../covariates.tsv"
bgen_file="../chr22.bgen"
sample_file="../chr22.sample"

regenie_gene_anno="annotation.tsv"
regenie_gene_setlist="set_list.tsv"
regenie_gene_masks="masks.tsv"

regenie \
    --step 2 \
    --bgen ${bgen_file} \
    --ref-first \
    --sample ${sample_file} \
    --phenoFile ${phenotype_file} \
    --bt --firth --approx --pThresh 0.01 \
    --covarFile ${covariates_file} \
    --bsize 400 \
    --pred fit_bin_l1_pred.list \
    --anno-file ${regenie_gene_anno} \
    --set-list ${regenie_gene_setlist} \
    --mask-def ${regenie_gene_masks} \
    --aaf-bins 0.01,0.001 \
    --threads 2 \
    --gz \
    --check-burden-files \
    --split \
    --out regenie_step2_out