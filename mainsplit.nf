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
  tuple val(phenotype), path(chunks)
  path master
  

  output:
  path("fit_bin_l1*"), emit: regenie_step1_l1_out
  //path("fit_bin_l1_${phenotype}_pred.list"), emit: regenie_step1_l1_predlist

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
    --out fit_bin_l1_${phenotype} \
    --l1-phenoList ${phenotype} \
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

    step1_l0_out = RegenieStep1_L0.out.regenie_step1_l0_out.flatten()

    // group by phenotype
    
    phenoMap_ch = phenotype_file.map { file ->
      def m = [:]
      file.withReader { r ->
          def header = r.readLine().split('\t')
          (2..<header.size()).each { i -> m["Y${i - 1}"] = header[i] }
      }
      return m
    }

    step1_l0_out_by_pheno = step1_l0_out.combine(phenoMap_ch)
        .map { file, phenoMap ->
            def group_key = file.getName().split('_')[-1]
            tuple(phenoMap[group_key], file)
        }
        .groupTuple().view()

/*
    // build map from Y_n to phenotype name
    def phenoMap = [:]
    def phenotypes_array = phenotype_file.get().newReader().readLine().split("\t")
    for (int i = 2; i < phenotypes_array.length; i++){
            phenoMap["Y" + i - 1] = phenotypes_array[i]
    }

    //println "phenoMap ${phenoMap}"

    step1_l0_out_by_pheno = step1_l0_out.collect().flatten()
      .map {file ->
            def name = file.getName()
            def group_key = name.split('_')[-1]
            tuple(phenoMap[group_key], file)
      }
      .groupTuple()
*/

    RegenieStep1_L1(genotypes_array_tuple, phenotype_file, covariates_file, step1_l0_out_by_pheno, step1_master)

    regenie_step1_out = RegenieStep1_L1.out.regenie_step1_l1_out.flatten().collect()

    emit:
    regenie_step1_out
}

process RegenieStep2 {
  
  tag "regenie_step2_${phenotype_file.baseName}"
  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  path step1_out_files
  path phenotype_file
  path covariates_file
  path bgen_file
  path sample_file
  path regenie_gene_anno
  path regenie_gene_setlist
  path regenie_gene_masks

  output:
  path("regenie_step2_out_${phenotype_file.baseName}_${bgen_file.baseName}_*.gz"), emit: regenie_step2_out

  script:
  def bt_flag     = params.phenotypes_binary_trait ? "--bt" : ""
  def firth_flag  = params.regenie_firth ? "--firth" : ""
  def approx_flag = params.regenie_firth_approx ? "--approx" : ""
  """
  cat fit_bin_l1_*_pred.list > fit_bin_l1_pred.list

  regenie \
    --step 2 \
    --bgen ${bgen_file} \
    --ref-first \
    --sample ${sample_file} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} ${firth_flag} ${approx_flag} --pThresh 0.01 \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step2} \
    --pred fit_bin_l1_pred.list \
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


process MergePerPhenotype {

  tag "merge_${phenofile.baseName}"

  publishDir "${params.outdir}", mode: 'move'

  input:
  path phenofile
  path all_result_files

  output:
  path("merged/*.regenie.gz")

  
  script:
  """
  mkdir merged

  # Extract phenotype names from header
  head -n 1 ${phenofile} | cut -f3- | tr '\t' '\n' > pheno_names.txt

  # Loop over phenotype names
  while read pheno; do
    # Extract header from the first matching file
    first_file=\$(ls ${all_result_files} | grep "_\${pheno}.regenie.gz" | head -n 1)
    zcat "\$first_file" | head -n 2 > "merged/\${pheno}.regenie"

    # Concatenate all matching files, skip their headers, and sort
    ls ${all_result_files} | grep "_\${pheno}.regenie.gz" | while read f; do
      zcat "\$f" | tail -n +3
    done | sort -k1,1 -k2,2n >> "merged/\${pheno}.regenie"

    # Compress
    bgzip -f "merged/\${pheno}.regenie"
  done < pheno_names.txt
  """
}



workflow Regenie {
  take:
  genotypes_array_tuple   // tuple val(plink_root), path(bed), path(bim), path(fam)
  phenotype_file          // path(phenotype_file)
  covariates_file         // path(covariates_file)
  bgen_files              // path(bgen_files)
  sample_file             // path(sample_file)

  main:
  RegenieStep1(genotypes_array_tuple, phenotype_file, covariates_file, 10)
  regenie_step1_out = RegenieStep1.out.regenie_step1_out.collect()

  regenie_anno_file    = file(params.regenie_gene_anno, checkIfExists: true)
  regenie_setlist_file = file(params.regenie_gene_setlist, checkIfExists: true)
  regenie_masks_file   = file(params.regenie_gene_masks, checkIfExists: true)

  bgen_file_ch = Channel.from(bgen_files).view()
  regenie_step1_out.view()
  RegenieStep2(regenie_step1_out, phenotype_file, covariates_file, bgen_file_ch, sample_file, 
               regenie_anno_file, regenie_setlist_file, regenie_masks_file)

  step2_out_files = RegenieStep2.out.regenie_step2_out.flatten().collect()
  MergePerPhenotype(phenotype_file, step2_out_files)

  emit:
  MergePerPhenotype.out
}

workflow {

  genotypes_array_ch = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
  genotypes_array_tuple = genotypes_array_ch.map{name, files -> tuple(name, files[1], files[0], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bin), path(fam)

  pheno_file_ch = Channel.fromPath(params.phenotypes_files).first()
  covariates_file = file(params.covariates_file)
  bgen_files = file(params.genotypes_bgen).findAll { !it.toString().contains("chrY") }
  sample_file = file(params.sample_file)

  Regenie(genotypes_array_tuple, pheno_file_ch, covariates_file, bgen_files, sample_file)

}
