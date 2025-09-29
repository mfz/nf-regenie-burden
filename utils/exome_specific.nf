
nextflow.enable.dsl=2

params.genome_bgens = "${WORKSPACE_BUCKET}/acaf_bgi"
params.exome_bgens = "${WORKSPACE_BUCKET}/exome_bgen"
params.exome_only_bgens = "${WORKSPACE_BUCKET}/exome_only"


process MakeExomeOnly {
    tag { chr }
    container "${ARTIFACT_REGISTRY_DOCKER_REPO}/florianzink/nf-gwas-gcloud:v0.3"
    scratch true
    disk '50GB'
    cpus 4
    memory '8 GB'
   	
    publishDir "${params.exome_only_bgens}"

    input:
    tuple val(chr), path(genome_bgi, stageAs: 'genome.bgi'), path(exome_bgen), path(exome_bgi)

    output:
    path("chr${chr}_exome_only_*.bgen")
    path("chr${chr}_exome_only_*.bgen.bgi")

    script:
    """
    # 1. extract variant_ids from sqlite indices
    sqlite3 ${genome_bgi} "SELECT rsid FROM Variant;" | sort -u > genome.ids
    sqlite3 ${exome_bgi}  "SELECT rsid FROM Variant;" | sort -u > exome.ids

    # 2. compute exome-only ids
    comm -23 exome.ids genome.ids > exome_only.ids

    # 3. split into chunks of 200,000 lines
    split -l 200000 -d --additional-suffix=.ids exome_only.ids chunk_

    # 4. loop over chunks and extract variants
    for ids in chunk_*.ids; do
        chunk=\$(basename \$ids .ids)
        out_bgen=chr${chr}_exome_only_\${chunk}.bgen
        bgenix -g ${exome_bgen} -incl-rsids \$ids > \$out_bgen
        bgenix -g \$out_bgen -index
    done
    """
}



workflow {
    Channel
        .fromPath("${params.genome_bgens}/chr*.bgen.bgi")
        .map { genome_bgi ->
            def chr = genome_bgi.name.replaceFirst(/chr([0-9XYMT]+).bgen.bgi/, "\$1")
            tuple(
                chr,
                file(genome_bgi),           // unique name
                file("${params.exome_bgens}/chr${chr}.bgen"),
                file("${params.exome_bgens}/chr${chr}.bgen.bgi")
            )
        }
        | MakeExomeOnly
}
