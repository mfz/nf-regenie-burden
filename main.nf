// Define Regenie subworkflow to encapsulate Step1, Step2, and merge per phenotype

nextflow.enable.dsl=2


process RegenieStep1 {

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  path phenotype_file   //value
  path covariates_file  // value

  output:
  path("regenie_step1_out*"), emit: regenie_step1_out

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
    --out regenie_step1_out
  """
}

process RegenieStep2 {

  input:
  path step1_out_files
  path phenotype_file
  path covariates_file
  path bgen_file
  path sample_file

  output:
  path("regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}_*.gz"), emit: regenie_step2_out

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
    --covarFile ${covariates_file} \
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
  genotypes_array_tuple   // tuple val(plink_root), path(bed), path(bim), path(fam)
  phenotype_file          // path(phenotype_file)
  covariates_file         // path(covariates_file)
  bgen_files              // path(bgen_files)
  sample_file             // path(sample_file)

  main:
  

  RegenieStep1(genotypes_array_tuple, phenotype_file, covariates_file)

  // path(regenie_step1_out*)
  step1_out_files = RegenieStep1.out.regenie_step1_out

  bgen_file_ch = Channel.from(bgen_files)
  RegenieStep2(step1_out_files, phenotype_file, covariates_file, bgen_file_ch, sample_file)

  // path(step2_out_files)
  // RegenieStep2 is called with a queue channel
  // need to flatten and collect
  step2_out_files = RegenieStep2.out.regenie_step2_out.flatten().collect()

  MergePerPhenotype(phenotype_file, step2_out_files)

  emit:
  MergePerPhenotype.out
}




workflow {

  genotypes_array_ch = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
  genotypes_array_tuple = genotypes_array_ch.map{name, files -> tuple(name, files[1], files[0], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bin), path(fam)

  pheno_file_ch = Channel.fromPath(params.phenotypes_files)
  covariates_file = file(params.covariates_file)
  bgen_files = Channel.fromPath(params.genotypes_bgen).filter { !it.toString().contains("chrY") }.collect()  // make this a list
  sample_file = file(params.sample_file)

  RegenieSubworkflow(genotypes_array_tuple, pheno_file_ch, covariates_file, bgen_files, sample_file)
}
