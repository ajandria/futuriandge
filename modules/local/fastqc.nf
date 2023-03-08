
process fastqc {

    publishDir "${params.outDir}/FastQC", mode: 'symlink'

    label = 'intense'

    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*"), emit: fastqc_on_raw_to_multiqc

    """
    fastqc -t ${task.cpus} ${reads}
    """  
}
