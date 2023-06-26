process SAMTOOLS_FLAGSTAT {
    tag "${meta.id}"

    publishDir "${params.outDir}/samtools", mode: 'symlink'

    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}_flagstat.txt"), emit: samtools_to_multiqc

    script:
    """
    samtools flagstat ${bam} -@ ${task.cpus} > ${meta.id}_flagstat.txt
    """
}
