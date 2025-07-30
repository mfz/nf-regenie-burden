nextflow.enable.dsl=2

params.bgen_dir = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/acaf_threshold/bgen/chr*{.bgen,.bgen.bgi}"  // Input BGEN directory in GCS
params.output_dir = "${WORKSPACE_BUCKET}/split_bgen" // Output directory for split files
params.chunk_size = 200000 

process split_bgen {
    cpus 1
    memory '2GB'
    disk '1000GB'

    publishDir "${params.output_dir}", mode: 'move'

    input:
    tuple val(bgen_filename), path(bgen_files)

    output:
    path "*_*_*.bgen"

    script:
    """

    sqlite3 -csv ${bgen_filename}.bgen.bgi "
        WITH variant_ordered AS (
        SELECT chromosome, position,
           ROW_NUMBER() OVER (PARTITION BY chromosome ORDER BY position) AS row_num
        FROM Variant
        ),
        chunked AS (
        SELECT chromosome,
           MIN(position) AS startpos,
           MAX(position) AS endpos
        FROM variant_ordered
        GROUP BY chromosome, (row_num - 1) / ${params.chunk_size}
        )
        SELECT chromosome, startpos, endpos FROM chunked ORDER BY chromosome, startpos;" > chunks.csv

    
    cat chunks.csv | while IFS=, read -r CHROM STARTPOS ENDPOS; do
    OUTPUT_FILE="\${CHROM}_\${STARTPOS}_\${ENDPOS}.bgen"
    
    # Run bgenix to extract the chunk
    bgenix -g "${bgen_filename}.bgen" -incl-range "\${CHROM}:\${STARTPOS}-\${ENDPOS}" > "\$OUTPUT_FILE"
    done
    """
} 

workflow {
    Channel.fromFilePairs("${params.bgen_dir}")
        | split_bgen
}
