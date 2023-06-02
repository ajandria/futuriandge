process picard_metrics {
    tag "${meta.id}"
    label 'process_high'

    publishDir "${params.outDir}/picard", mode: 'symlink'
   
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}_picard_CollectRnaSeqMetrics"), emit: picard_to_multiqc

    script:
    // Dict of library types
    def strandedness_opts = ['reverse', 'forward', 'unstranded']
    
    // By default, try to auto-detect library type
    def strandedness = 'NONE'

    // Check if library type is valid and assign library type
    if (strandedness_opts.contains(meta.strandedness)) {
        if (meta.strandedness == 'reverse') {
            strandedness = 'SECOND_READ_TRANSCRIPTION_STRAND'
        } else if (meta.strandedness == 'forward') {
            strandedness = 'FIRST_READ_TRANSCRIPTION_STRAND'
        } else if (meta.strandedness == 'unstranded') {
            strandedness = 'NONE'
        }
    } else {
        log.info "[picard CollectRnaSeqMetrics] Invalid strandedness specified '--libType=${lib_type}', defaulting to 'NONE'."
    }

    """
    picard CollectRnaSeqMetrics \
        -I ${bam} \
        -O ${meta.id}_picard_CollectRnaSeqMetrics \
        --REF_FLAT ${params.refFlat} \
        --STRAND ${strandedness}
    """
}
