
process fastp {

    publishDir "${params.outDir}/fastp", mode: 'symlink'

    label = 'intense'

    tag "fastp on ${meta}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_fastp_R*"), emit: fastp_to_align
    path("*"), emit: fastp_to_multiqc

    script:
        """
        cp ${reads[0]} ${meta.id}_raw_R1.fastq.gz
        cp ${reads[1]} ${meta.id}_raw_R2.fastq.gz
        fastp -i ${meta.id}_raw_R1.fastq.gz -I ${meta.id}_raw_R2.fastq.gz \
            -o ${meta.id}_fastp_R1.fastq.gz -O ${meta.id}_fastp_R2.fastq.gz \
            -j ${meta.id}_fastp.json \
            -h ${meta.id}_fastp.html \
            -q 30 \
            -e 25 \
            -n 5 \
            -l 40 \
            -c \
            -x \
            -p \
            --verbose \
            -w ${task.cpus}
        rm ${meta.id}_raw_R1.fastq.gz
        rm ${meta.id}_raw_R2.fastq.gz
        """
}