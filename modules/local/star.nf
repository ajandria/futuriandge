process star {
    publishDir "${params.outDir}/star", mode: 'symlink'
    label = 'intense'
    tag "star on ${meta}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' :
        'biocontainers/mulled-v2-1fa26d1ce03c295fe2fdcf85831a92fbcbd7e8c2:1df389393721fc66f3fd8778ad938ac711951107-0' }"


    input:
    tuple val(meta), path(reads)
    
    output:
    tuple val(meta), path("${meta.id}*.bam"), emit: aligned_with_star
    path("*"), emit: star_to_multiqc
    
    script:
    def input_reads = meta.single_end == 'true' ? "--readFilesIn ${reads[0]}" : "--readFilesIn ${reads[0]} ${reads[1]}"
    
    """
    STAR --runMode alignReads \\
    --genomeDir ${params.reference_genome} \\
    --sjdbGTFfile ${params.gtf} \\
    ${input_reads} \\
    --readFilesCommand zcat \\
    --outSAMtype BAM SortedByCoordinate \\
    --outFileNamePrefix ${meta.id} \\
    --runThreadN ${task.cpus}
    
    samtools index -@ ${task.cpus} *.bam
    """
}
