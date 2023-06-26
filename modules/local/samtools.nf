process samtools {
    tag "${meta.id}"

    publishDir "${params.outDir}/samtools", mode: 'symlink'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0' :
        'biocontainers/samtools:1.17--h00cdaf9_0' }"


    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}_flagstat.txt"), emit: samtools_to_multiqc

    script:
    """
    samtools flagstat ${bam} -@ ${task.cpus} > ${meta.id}_flagstat.txt
    """
}
