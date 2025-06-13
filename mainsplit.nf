nextflow.enable.dsl=2


process RegenieStep1_Split {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  tuple val(meta), path(phenotype_file)  
  path covariates_file  // value
  val num_chunks

  
  output:
  tuple val(meta), path("fit_bin_parallel*.snplist"), path("fit_bin_parallel.master"), emit: regenie_step1_split_out
  
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
  val meta
  path phenotype_file   //value
  path covariates_file  // value
  path step1_snplists
  path step1_master
  val job

  output:
  tuple val(meta), path("fit_bin_parallel*"), emit: regenie_step1_l0_out

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
    --threads ${task.cpus} \
    --out fit_bin_l0_${job} \
    --run-l0 fit_bin_parallel.master,${job}

  cp fit_bin_parallel.master copy.master
  """
}




process RegenieStep1_L1 {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  val meta
  path phenotype_file   //value
  path covariates_file  // value
  path locos
  path master
  val phenonum
  

  output:
  tuple val(meta), path("fit_bin_l1*"), emit: regenie_step1_l1_out
  //path("fit_bin_l1_${phenotype}_pred.list"), emit: regenie_step1_l1_predlist

  script:
  def bt_flag = params.phenotypes_binary_trait ? "--bt" : ""
  """
  phenotype=\$(head -n 1 ${phenotype_file} | cut -f\$((${phenonum}+2)) )

  regenie \
    --step 1 \
    --bed ${genotype_array} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step1} \
    --ref-first \
    --lowmem \
    --out fit_bin_l1_\${phenotype} \
    --l1-phenoList \${phenotype} \
    --run-l1 ${master} \
    --threads ${task.cpus} \
    --use-relative-path
  """
}


process RegenieStep2_Burden {
  
  tag "regenie_step2_${phenotype_file.baseName}"
  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  path step1_out_files
  val meta
  path phenotype_file
  path covariates_file
  path bgen_file
  path sample_file
  path regenie_gene_anno
  path regenie_gene_setlist
  path regenie_gene_masks

  output:
  tuple val(meta), path("regenie_step2_out_${bgen_file.baseName}_*.gz"), emit: regenie_step2_out

  script:
  def bt_flag     = params.phenotypes_binary_trait ? "--bt" : ""
  def firth_flag  = params.regenie_firth ? "--firth --firth-se --pThresh 0.05" : ""
  def approx_flag = params.regenie_firth_approx ? "--approx" : ""
  """
  cat fit_bin_l1_*_pred.list > fit_bin_l1_pred.list

  regenie \
    --step 2 \
    --bgen ${bgen_file} \
    --ref-first \
    --sample ${sample_file} \
    --phenoFile ${phenotype_file} \
    ${bt_flag} ${firth_flag} ${approx_flag}  \
    --covarFile ${covariates_file} \
    --bsize ${params.regenie_bsize_step2} \
    --pred fit_bin_l1_pred.list \
    --anno-file ${regenie_gene_anno} \
    --set-list ${regenie_gene_setlist} \
    --mask-def ${regenie_gene_masks} \
    --aaf-bins ${params.regenie_gene_aaf} \
    --threads ${task.cpus} \
    --gz \
    --check-burden-files \
    --out regenie_step2_out_${bgen_file.baseName}
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




workflow {

  genotypes_array_tuple = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
                          .map{name, files -> tuple(name, files[1], files[0], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bin), path(fam)

  pheno_file_ch = Channel.fromPath(params.phenotypes_files) 
                  .map {file ->
                        def meta = file.baseName
                        [meta, file]
                  }
  // channel of tuple val(meta), path(pheno_file)

  covariates_file = file(params.covariates_file)

  RegenieStep1_Split(genotypes_array_tuple, 
                     pheno_file_ch, 
                     covariates_file, 
                     10)

  step1_split_out = RegenieStep1_Split.out.regenie_step1_split_out
  // channel tuple val(meta), path(*.snplist), path(master)

  step1_l0_in = pheno_file_ch  // val(meta), path(pheno_file)
      .join(step1_split_out)   // val(meta), path(*.snplist), path(master)


  // scatter into 10 jobs
  jobs = Channel.from(1..10)

  combined_step1_l0_in = step1_l0_in.combine(jobs)

  RegenieStep1_L0(genotypes_array_tuple,
                                 combined_step1_l0_in.map {it[0]}, // meta
                                 combined_step1_l0_in.map {it[1]}, // phenofile
                                 covariates_file,
                                 combined_step1_l0_in.map {it[2]}, // snplists
                                 combined_step1_l0_in.map {it[3]}, // master
                                 combined_step1_l0_in.map {it[4]}) // job
  
  step1_l0_out = RegenieStep1_L0.out.regenie_step1_l0_out
  // channel tuple val(meta), path(jobJ_Y*)

  // gather
  // need to group by val(meta)
  step1_l0_out_grouped = step1_l0_out
    .groupTuple()
    .map {key, val ->
          def flat = val.flatten()
          tuple(key, flat)}
  // channel val(meta), path(job*_Y*)

  // scatter by phenotype
  // input to RegenieStep1_L1 is
  // geno, pheno , covar, job*_Y*, master, phenonum
  step1_l1_in = pheno_file_ch      // val(meta), path(pheno)
      .join(step1_split_out)       // val(meta), path(*snplist), path(master)
      .join(step1_l0_out_grouped)  // val(meta), path(job*_Y*)

  phenonums = Channel.from(1..params.num_phenotypes_per_file)
  combined_step1_l1_in = step1_l1_in.combine(phenonums)

  RegenieStep1_L1(genotypes_array_tuple,
                                 combined_step1_l1_in.map {it[0]}, // meta
                                 combined_step1_l1_in.map {it[1]}, // phenotype_file
                                 covariates_file,
                                 combined_step1_l1_in.map {it[4]}, // job*_Y*
                                 combined_step1_l1_in.map {it[3]}, // master
                                 combined_step1_l1_in.map {it[5]}) // phenonum
  
  step1_l1_out = RegenieStep1_L1.out.regenie_step1_l1_out
  
  // channel tuple val(meta), path(fit_bin_l1_pheno)

  // gather phenotypes again
  step1_l1_out_grouped = step1_l1_out
    .groupTuple()
    .map {key, val ->
          def flat = val.flatten()
          tuple(key, flat)}
  // channel tuple val(meta), path(fit_bin_l1_*)

  bgen_ch = Channel.fromPath(params.genotypes_bgen)
     .filter {file -> !file.toString().contains("chrY")}

  sample_file = file(params.sample_file)
  regenie_anno_file    = file(params.regenie_gene_anno, checkIfExists: true)
  regenie_setlist_file = file(params.regenie_gene_setlist, checkIfExists: true)
  regenie_masks_file   = file(params.regenie_gene_masks, checkIfExists: true)

  // scatter over bgens
  combined_step2_in = pheno_file_ch           // val(meta), path(pheno_file)        
      .join(step1_l1_out_grouped)             // val(meta), path(fit_bin_l1_*)
      .combine(bgen_ch)                       // path(bgen_file)


  RegenieStep2_Burden(combined_step2_in.map {it[2]}, // path(fit_bin_l1_*)
                combined_step2_in.map {it[0]}, // val(meta)
                combined_step2_in.map {it[1]}, // path(pheno_file)
                covariates_file,
                combined_step2_in.map {it[3]}, // path(bgen_file)
                sample_file,
                regenie_anno_file,
                regenie_setlist_file,
                regenie_masks_file)

  step2_out = RegenieStep2_Burden.out.regenie_step2_out

  // gather over bgens
  step2_out_grouped = step2_out
    .groupTuple()
    .map {key, val ->
          def flat = val.flatten()
          tuple(key, flat)}
  
  merge_in = pheno_file_ch     // val(meta), path(pheno_file)
    .join(step2_out_grouped)   // val(meta), path(step2out)

  MergePerPhenotype( merge_in.map {it[1]},
                     merge_in.map {it[2]})

}
