nextflow.enable.dsl=2


process RegenieStep1_Split {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  path phenotype_file   //value
  path covariates_file  // value
  val num_chunks

  
  output:
  path("fit_bin_parallel*"), emit: regenie_step1_split_out
  path("fit_bin_parallel.master"), emit: regenie_step1_master

  script:
  def bt_flag = params.phenotypes_binary_trait ? "--bt" : ""
  """
  regenie \
    --step 1 \
    --bed ${genotype_array} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step1} \
    --ref-first \
    --lowmem \
    --out fit_bin_l0 \
    --split-l0 fit_bin_parallel,${num_chunks}
  """
}


process RegenieStep1_L0 {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  path phenotype_file   //value
  path covariates_file  // value
  path step1_split
  val job

  output:
  path("fit_bin_parallel*"), emit: regenie_step1_l0_out

  script:
  def bt_flag = params.phenotypes_binary_trait ? "--bt" : ""
  """
  regenie \
    --step 1 \
    --bed ${genotype_array} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step1} \
    --ref-first \
    --lowmem \
    --out fit_bin_l0_${job} \
    --run-l0 fit_bin_parallel.master,${job}
  """
}


process RegenieStep1_L1 {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  path phenotype_file   //value
  path covariates_file  // value
  path step1_l0
  path master
  

  output:
  path("fit_bin_l1*"), emit: regenie_step1_l1_out

  script:
  def bt_flag = params.phenotypes_binary_trait ? "--bt" : ""
  """
  regenie \
    --step 1 \
    --bed ${genotype_array} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step1} \
    --ref-first \
    --lowmem \
    --out fit_bin_l1 \
    --run-l1 ${master} \
    --use-relative-path
  """
}


workflow RegenieStep1 {

    take:
    genotypes_array_tuple   // tuple val(plink_root), path(bed), path(bim), path(fam)
    phenotype_file          // path(phenotype_file)
    covariates_file         // path(covariates_file)
    num_chunks

    main:
    RegenieStep1_Split(genotypes_array_tuple, phenotype_file, covariates_file, num_chunks)
    step1_split_out = RegenieStep1_Split.out.regenie_step1_split_out.collect()
    step1_master = RegenieStep1_Split.out.regenie_step1_master

    jobs = Channel.from(1..num_chunks)
    RegenieStep1_L0(genotypes_array_tuple, phenotype_file, covariates_file, step1_split_out, jobs)

    step1_l0_out = RegenieStep1_L0.out.regenie_step1_l0_out.collect()
    RegenieStep1_L1(genotypes_array_tuple, phenotype_file, covariates_file, step1_l0_out, step1_master)

    regenie_step1_out = RegenieStep1_L1.out.regenie_step1_l1_out

    emit:
    regenie_step1_out
}


workflow {

  genotypes_array_ch = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
  genotypes_array_tuple = genotypes_array_ch.map{name, files -> tuple(name, files[1], files[0], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bin), path(fam)

  pheno_file_ch = Channel.fromPath(params.phenotypes_files).first()
  covariates_file = file(params.covariates_file)

  RegenieStep1(genotypes_array_tuple, pheno_file_ch, covariates_file, 10)

}
