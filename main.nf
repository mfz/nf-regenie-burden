// Define Regenie subworkflow to encapsulate Step1, Step2, and merge per phenotype

nextflow.enable.dsl=2



process RegenieStep1 {

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam) 
  tuple val(meta), path(phenotype_file)   
  path covariates_file  

  output:
  tuple val(meta), path(phenotype_file), path("regenie_step1_*"), emit: regenie_step1_out

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
    --out regenie_step1_${meta.phenotype}
  """
}

process RegenieStep2 {

  input:
  val(meta)
  path phenotype_file
  path step1_out_files
  path covariates_file
  path bgen_file
  path sample_file
  path regenie_gene_anno
  path regenie_gene_setlist
  path regenie_gene_masks

  output:
  tuple val(meta), path("regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}_*.gz"), emit: regenie_step2_out

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
    --anno-file ${regenie_gene_anno} \
    --set-list ${regenie_gene_setlist} \
    --mask-def ${regenie_gene_masks} \
    --aaf-bins ${params.regenie_gene_aaf} \
    --threads 2 \
    --gz \
    --check-burden-files \
    --out regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}
  """
}


process MergePhenotype {


  publishDir "${params.outdir}", mode: 'move'

  input:
  tuple val(phenotype), path(step2_out_files)
  

  output:
  path("merged/*.regenie.gz")

  
  script:
  """
  mkdir merged

  # Extract header from the first matching file
  first_file=\$(ls | grep ".regenie.gz" | head -n 1)
  zcat "\$first_file" | head -n 2 > "merged/${phenotype}.regenie"

  # Concatenate all matching files, skip their headers, and sort
  ls | grep ".regenie.gz" | while read f; do
      zcat "\$f" | tail -n +3
    done | sort -k1,1 -k2,2n >> "merged/${phenotype}.regenie"

  # Compress
  bgzip -f "merged/${phenotype}.regenie"
  """
}




workflow {

  genotypes_array_tuple = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
                          .map{name, files -> tuple(name, files[1], files[0], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bin), path(fam)

  pheno_file_ch = Channel.fromPath(params.phenotypes_files) 
                  .map {file ->
                        def meta = [phenotype: file.baseName]
                        [meta, file]
                  }.view()
  // channel of tuple val(meta), path(pheno_file)

  covariates_file = file(params.covariates_file)

  RegenieStep1(genotypes_array_tuple, pheno_file_ch, covariates_file)
  step1_out_ch = RegenieStep1.out.regenie_step1_out  // tuple val(meta), path(pheno_file), path(regenie_step1_*)

  bgen_files_ch = Channel
      .fromPath(params.genotypes_bgen, checkIfExists: true)
      .filter { !it.toString().contains("chrY") }
  
  sample_file = file(params.sample_file)

  combined_ch = step1_out_ch.combine(bgen_files_ch).view()

  regenie_anno_file    = file(params.regenie_gene_anno, checkIfExists: true)
  regenie_setlist_file = file(params.regenie_gene_setlist, checkIfExists: true)
  regenie_masks_file   = file(params.regenie_gene_masks, checkIfExists: true)

  RegenieStep2(combined_ch.map {it[0]},
               combined_ch.map {it[1]},
               combined_ch.map {it[2]},
               covariates_file,
               combined_ch.map {it[3]},
               sample_file,
               regenie_anno_file,
               regenie_setlist_file,
               regenie_masks_file)

  step2_out_ch = RegenieStep2.out.regenie_step2_out  // tuple val(meta), path(regenie_step2_out*.gz)

  // group by phenotyoe
  grouped_ch = step2_out_ch
      .flatMap { meta, files -> 
                 files.collect { file -> tuple(meta.phenotype, file) }
      }
  .groupTuple().view()

  MergePhenotype(grouped_ch)


}
