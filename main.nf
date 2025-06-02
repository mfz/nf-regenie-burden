// Define Regenie subworkflow to encapsulate Step1, Step2, and merge per phenotype

nextflow.enable.dsl=2


process RegenieStep1 {

  input:
  path genotype_array
  path phenotype_file
  path covariates_file

  output:
  tuple path(phenotype_file), path("regenie_step1_out*"), emit: regenie_step1_out

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
    --lowmem \
    --out regenie_step1_out
  """
}

process RegenieStep2 {

  

  input:
  tuple path(phenotype_file), path(step1_out)
  path bgen_file
  path sample_file

  output:
  tuple path(phenotype_file), path("regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}_*.gz"), emit: regenie_step2_out

  script:
  def bt_flag     = params.phenotypes_binary_trait ? "--bt" : ""
  def firth_flag  = params.regenie_firth ? "--firth" : ""
  def approx_flag = params.regenie_firth_approx ? "--approx" : ""
  """
  regenie \
    --step 2 \
    --bgen ${bgen_file} \
    --ref-first \
    --sample ${sample_file} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} ${firth_flag} ${approx_flag} --pThresh 0.01 \
    --bsize ${params.regenie_bsize_step2} \
    --pred regenie_step1_out_pred.list \
    --anno-file ${params.regenie_gene_anno} \
    --set-list ${params.regenie_gene_setlist} \
    --mask-def ${params.regenie_gene_masks} \
    --aaf-bins ${params.regenie_gene_aaf} \
    --threads 2 \
    --gz \
    --check-burden-files \
    --split \
    --out regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}
  """
}

process MergePerPhenotype {

  cpus 2
  memory '4G'

  tag "merge_${phenofile.baseName}"

  publishDir "${params.outdir}/results", mode: 'move'

  input:
  path phenofile
  path all_result_files

  output:
  path("*.txt.gz")

  
  script:
  """
  head -n 1 ${phenofile} | cut -f3- | tr '\t' '\n' > pheno_names.txt
  
  while read pheno; do
    zcat \$(ls ${all_result_files} | grep "_\${pheno}.gz") | sort -k1,1 -k2,2n | bgzip -c > \${pheno}.txt.gz
  done < pheno_names.txt
  """
}

workflow RegenieSubworkflow {
  take:
  phenotype_file
  bgen_files

  main:
  def genotype_array = file(params.genotypes_array)
  def covariates_file = file(params.covariates_file)
  def sample_file = file(params.sample_file)


  RegenieStep1(genotype_array, phenotype_file, covariates_file)

  step1_out = RegenieStep1.out.regenie_step1_out.first()

  RegenieStep2(step1_out, bgen_files, sample_file)

  step2_out = RegenieStep2.out.regenie_step2_out.collect()

  MergePerPhenotype(phenotype_file, step2_out)

  emit:
  MergePerPhenotype.out
}




workflow {
  Channel.fromPath(params.phenotypes_files).set { pheno_files }
  Channel.fromPath(params.genotypes_bgen).set { bgen_files }

  RegenieSubworkflow(pheno_files.first(), bgen_files)
}
