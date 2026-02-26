nextflow.enable.dsl=2


process PlinkMacFilter {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)
  tuple val(meta), path(phenotype_file)

  output:
  tuple val(meta), path("${meta}.snplist"), emit: plink_mac_snplist

  script:
  def mac_threshold = params.plink_mac ?: 5
  def maf_threshold = params.plink_maf ?: 0.01
  """
  # Extract samples from phenotype file (FID IID only)
  cut -f1,2 ${phenotype_file} | tail -n +2 > samples.txt

  # Run PLINK with MAC filter and keep only these samples
  plink \
    --bfile ${genotype_array} \
    --keep samples.txt \
    --mac ${mac_threshold} \
    --write-snplist \
    --out ${meta}
  """
}



process RegenieStep1_Split {

  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(genotype_array), path(plink_bed), path(plink_bim), path(plink_fam)   // value
  tuple val(meta), path(phenotype_file), path(mac_snplists)  // value  
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
    --extract ${mac_snplists} \
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
  path covariates_file // value
  tuple val(meta), path(phenotype_file), path(mac_snplists), path(step1_snplists), path(step1_master), val(job) // channel
  
  output:
  tuple val(meta), path("fit_bin_parallel*"), emit: regenie_step1_l0_out

  script:
  def bt_flag = params.phenotypes_binary_trait ? "--bt" : ""
  """
  regenie \
    --step 1 \
    --bed ${genotype_array} \
    --extract ${mac_snplists} \
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
  path covariates_file // value
  tuple val(meta), path(phenotype_file), path(mac_snplists), path(step1_snplists), path(master), path(locos), val(phenonum) // channel
 

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
    --extract ${mac_snplists} \
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
  tuple val(meta), path(phenotype_file), path(mac_snplists), path(step1_out_files), path(bgen_file)
  path covariates_file
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


process RegenieStep2 {
  
  tag "regenie_step2_${phenotype_file.baseName}"
  publishDir "${params.outdir}/logs", pattern: "*.log", mode: "copy"

  input:
  tuple val(meta), path(phenotype_file), path(mac_snplists), path(step1_out_files), path(bgen_file)
  path covariates_file
  path sample_file

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
    --threads ${task.cpus} \
    --gz \
    --out regenie_step2_out_${bgen_file.baseName}
  """
}





process MergePerPhenotype {

  tag "merge_${pheno}"

  publishDir "${params.outdir}", mode: 'move'

  input:
  val pheno
  path result_files

  output:
  path("merged/*.regenie.gz")


  script:
  """
  mkdir merged

  # Extract header from the first matching file ignoring lines starting with #
  first_file=\$(ls ${result_files} | head -n 1)
  zcat "\$first_file" | head -n 2 | grep -v '^#' | head -n 1 > "merged/${pheno}.regenie"
  
    # Concatenate all matching files, skip lines starting with # and one header line, and sort
  ls ${result_files} | while read f; do
      zcat "\$f" | grep -v '^#' | tail -n +2
  done | sort -k1,1 -k2,2n >> "merged/${pheno}.regenie"


  # Compress
  bgzip -f "merged/${pheno}.regenie"

  """
}


workflow {

  genotypes_array_tuple = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
                          .map{name, files -> tuple(name, files[0], files[1], files[2])}.first()
  // tuple val(plink_root), path(bed), path(bim), path(fam)  

  pheno_file_ch = Channel.fromPath(params.phenotypes_files) 
                  .map {file ->
                        def meta = file.baseName
                        [meta, file]
                  }
  // channel of tuple val(meta), path(pheno_file)

  covariates_file = file(params.covariates_file)


  PlinkMacFilter(genotypes_array_tuple,
                 pheno_file_ch)

  mac_snplists = PlinkMacFilter.out.plink_mac_snplist
  // tuple val(meta), path(${meta}.snplist)

  pheno_with_snplist = pheno_file_ch
  .join(mac_snplists)  // val(meta), path(pheno_file), path(snplist)



  RegenieStep1_Split(genotypes_array_tuple, 
                     pheno_with_snplist, 
                     covariates_file, 
                     10)

  step1_split_out = RegenieStep1_Split.out.regenie_step1_split_out
  // channel tuple val(meta), path(*.snplist), path(master)

  step1_l0_in = pheno_with_snplist // val(meta), path(pheno_file), path(macsnplist)
      .join(step1_split_out)   // val(meta), path(*.snplist), path(master)


  // scatter into 10 jobs
  jobs = Channel.from(1..10)

  combined_step1_l0_in = step1_l0_in.combine(jobs) // val(meta), path(pheno_file), path(macsnplist), path(*.snplist), path(master), val(job)

  RegenieStep1_L0(genotypes_array_tuple,
                 covariates_file,
                 combined_step1_l0_in)
  
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
  step1_l1_in = pheno_with_snplist      // val(meta), path(pheno), path(macsnplist)
      .join(step1_split_out)       // val(meta), path(*snplist), path(master)
      .join(step1_l0_out_grouped)  // val(meta), path(job*_Y*)

  phenonums = Channel.from(1..params.num_phenotypes_per_file)
  combined_step1_l1_in = step1_l1_in.combine(phenonums)
  // val(meta), path(pheno), path(macsnplist), path(*snplist), path(master), path(job*_Y*), val(phenonum)

  RegenieStep1_L1(genotypes_array_tuple,
                  covariates_file,
                  combined_step1_l1_in)



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

  // scatter over bgens
  combined_step2_in = pheno_with_snplist           // val(meta), path(pheno_file), path(mac_snplist)        
      .join(step1_l1_out_grouped)             // val(meta), path(fit_bin_l1_*)
      .combine(bgen_ch)                       // path(bgen_file)

  if (params.burden) {

   regenie_anno_file    = file(params.regenie_gene_anno, checkIfExists: true)
   regenie_setlist_file = file(params.regenie_gene_setlist, checkIfExists: true)
   regenie_masks_file   = file(params.regenie_gene_masks, checkIfExists: true)


    RegenieStep2_Burden(combined_step2_in
                covariates_file,
                sample_file,
                regenie_anno_file,
                regenie_setlist_file,
                regenie_masks_file)

    step2_out = RegenieStep2_Burden.out.regenie_step2_out
    // val(meta), path(step2_out_bgen_trait*)
  } else {

    RegenieStep2(combined_step2_in,
                covariates_file,
                sample_file)

    step2_out = RegenieStep2.out.regenie_step2_out
    // val(meta), path(step2_out_bgen_trait*)
  }
  // gather over bgens
  step2_out_grouped = step2_out
    .groupTuple()
    .map {key, val ->
          def flat = val.flatten()
          tuple(key, flat)}
  
  // step2_out_grouped is
  // val(meta), path(step2_out_bgen_trait*)

  merge_in = pheno_file_ch     // val(meta), path(pheno_file)
    .join(step2_out_grouped)   // val(meta), path(step2out)
    .flatMap {meta_, pheno_file_, step2out_ ->
        def header = pheno_file_.head(1).readLines().first().split(/\s+/)
        def phenos = header - ['FID','IID']
 
        phenos.collect { pheno ->
            // collect ALL matching files for this phenotype
            def matches = step2out_.findAll { it.name.endsWith("_${pheno}.regenie.gz") }
            matches ? tuple(meta_, pheno_file_, pheno, matches) : null
        }.findAll()
    }
    
MergePerPhenotype (merge_in.map {it[2]},
                   merge_in.map {it[3]})

}
