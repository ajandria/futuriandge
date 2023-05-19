
process star {

    publishDir "${params.outDir}/star", mode: 'symlink'

    label = 'intense'

    tag "star on ${meta}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}*.bam"), emit: aligned_with_star
    path("*"), emit: star_to_multiqc

    script:
     // Define input reads
     def input_reads = meta.single_end == 'true' ? "--readFilesIn ${reads[0]}" : "--readFilesIn ${reads[0]} ${reads[1]}"
        """
        STAR --runMode alignReads \
        --genomeDir ${params.reference_genome} \
        --sjdbGTFfile ${params.gtf} \
        --readFilesIn ${input_reads} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outFileNamePrefix ${meta.id} \
        --runThreadN ${task.cpus}

        samtools index -@ ${task.cpus} *.bam
        """  
}
