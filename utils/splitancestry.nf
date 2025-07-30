
params.baseDir = "${WORKSPACE_BUCKET}"
params.outDir = "${WORKSPACE_BUCKET}"
params.genotypes_array = "gs://fc-aou-datasets-controlled/v8/microarray/plink/arrays.{bed,bim,fam}"

// tab-delimited file with columns FID, IID, ancestry
params.ancestry_file = "${params.baseDir}/ancestry.tsv"
params.ancestries = ['AFR','AMR', 'EAS']
params.gender_file = "${params.baseDir}/is_male.tsv"

//SNP_PRUNING process
params.prune_enabled                         = true
params.prune_maf                             = 0.01
params.prune_window_kbsize                   = 1000
params.prune_step_size                       = 100
params.prune_r2_threshold                    = 0.9

//QC_FILTER process
params.qc_maf                                = 0.01
params.qc_mac                                = 100
params.qc_geno                               = 0.1
params.qc_hwe                                = '1e-15'
params.qc_mind                               = 0.1




process qc_and_filter_array {

    container  "${ARTIFACT_REGISTRY_DOCKER_REPO}/florianzink/nf-gwas-gcloud:v0.2"
    publishDir "${params.outDir}/arrays/${ancestry}/", mode: 'move'

    input:
    tuple val(array_filename), path(array_file)
    path(ancestry_file)
    path(gender_file)
    each ancestry

    output:
    path "${array_filename}_${ancestry}_qcfiltered.*"
    path "${array_filename}_${ancestry}_male_qcfiltered.*"
    path "${array_filename}_${ancestry}_female_qcfiltered.*"

    script:
    def prune = params.prune_enabled ? "--indep-pairwise ${params.prune_window_kbsize} ${params.prune_step_size} ${params.prune_r2_threshold}" : ''
    def extract = params.prune_enabled ? '--extract qcfiltered.prune.in' : '--extract qcfiltered.snplist'
    """
    awk '\$3==tolower("${ancestry}") {print 0,\$2}' "${ancestry_file}" > keep.ids

     plink2 \
     --bfile "${array_filename}" \
     --not-chr Y \
     --keep keep.ids \
     --maf ${params.qc_maf} \
     --mac ${params.qc_mac} \
     --geno ${params.qc_geno} \
     --hwe ${params.qc_hwe} \
     --mind ${params.qc_mind} \
     --write-snplist --write-samples --no-id-header \
     --out qcfiltered \
    $prune

    plink2 \
    --bfile "${array_filename}" \
    --keep qcfiltered.id \
    $extract \
    --make-bed \
    --out "${array_filename}_${ancestry}_qcfiltered"

    awk '{\$1=\$2; print \$0}' \
      "${array_filename}_${ancestry}_qcfiltered.fam" \
    > "${array_filename}_${ancestry}_qcfiltered.fam_fixed"

    mv "${array_filename}_${ancestry}_qcfiltered.fam_fixed" \
      "${array_filename}_${ancestry}_qcfiltered.fam"


    awk '\$3==1 {print \$1,\$2}' "${gender_file}" > male.ids
    awk '\$3==0 {print \$1,\$2}' "${gender_file}" > female.ids

    plink2 \
    --bfile "${array_filename}_${ancestry}_qcfiltered" \
    --keep male.ids \
    --maf ${params.qc_maf} \
    --mac ${params.qc_mac} \
    --geno ${params.qc_geno} \
    --hwe ${params.qc_hwe} \
    --mind ${params.qc_mind} \
    --make-bed \
    --out "${array_filename}_${ancestry}_male_qcfiltered"


    plink2 \
    --bfile "${array_filename}_${ancestry}_qcfiltered" \
    --keep female.ids \
    --maf ${params.qc_maf} \
    --mac ${params.qc_mac} \
    --geno ${params.qc_geno} \
    --hwe ${params.qc_hwe} \
    --mind ${params.qc_mind} \
    --make-bed \
    --out "${array_filename}_${ancestry}_female_qcfiltered"


    """
}

workflow {

    genotypes_array_ch = Channel.fromFilePairs(params.genotypes_array, size: 3, checkIfExists: true)
    ancestry_file_ch = Channel.fromPath(params.ancestry_file, checkIfExists: true)
    gender_file_ch = Channel.fromPath(params.gender_file, checkIfExists: true)    
    ancestry_ch = Channel.from(params.ancestries)

    qc_and_filter_array(
        genotypes_array_ch,
        ancestry_file_ch,
        gender_file_ch,
        ancestry_ch
    )

}

