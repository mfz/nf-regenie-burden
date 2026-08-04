
WORKSPACE_BUCKET="gs://workspace-bucket-wb-radiant-cabbage-3726"

nextflow.enable.dsl=2

params.genome_bgens = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/acaf_threshold/bgen"
params.exome_bgens = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/exome/bgen"
params.exome_only_bgens = "${WORKSPACE_BUCKET}/exome_only"


process MakeExomeOnly {
    tag { chr }
    container "florianzink/nf-gwas-gcloud:v0.3"
    scratch true
    disk '50GB'
    cpus 4
    memory '8 GB'

    publishDir "${params.exome_only_bgens}"

    input:
    tuple val(chr), path(genome_bgi, stageAs: 'genome.bgi'), path(exome_bgen), path(exome_bgi)

    output:
    path("chr${chr}_exome_only.bgen")
    path("chr${chr}_exome_only.bgen.bgi")

    script:
    """
    # 1. extract variant_ids from sqlite indices
    sqlite3 ${genome_bgi} "SELECT rsid FROM Variant;" | sort -u > genome.ids
    sqlite3 ${exome_bgi}  "SELECT rsid FROM Variant;" | sort -u > exome.ids

    # 2. compute exome-only ids
    comm -23 exome.ids genome.ids > exome_only.ids

    # 4. loop over chunks and extract variants
    out_bgen=chr${chr}_exome_only.bgen
    bgenix -g ${exome_bgen} -incl-rsids exome_only.ids > \$out_bgen
    bgenix -g \$out_bgen -index
    """
}



workflow {
    Channel
        .fromPath("${params.genome_bgens}/acaf_threshold.chr*.bgen.bgi")
        .map { genome_bgi ->
            def chr = genome_bgi.name.replaceFirst(/acaf_threshold.chr([0-9XYMT]+).bgen.bgi/, "\$1")
            tuple(
                chr,
                file(genome_bgi),           // unique name
                file("${params.exome_bgens}/exome.chr${chr}.bgen"),
                file("${params.exome_bgens}/exome.chr${chr}.bgen.bgi")
            )
        }
        | MakeExomeOnly
}

