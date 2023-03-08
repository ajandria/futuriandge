
process fastp {

    publishDir "${params.outDir}/fastp", mode: 'symlink'

    label = 'intense'

    tag "fastp on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_fastp_R*"), emit: fastp_to_align
    path("*"), emit: fastp_to_multiqc

    """
    cp ${reads[0]} ${sample_id}_raw_R1.fastq.gz
    cp ${reads[1]} ${sample_id}_raw_R2.fastq.gz
    fastp -i ${sample_id}_R1.fastq.gz -I ${sample_id}_R2.fastq.gz \
        -o ${sample_id}_fastp_R1.fastq.gz -O ${sample_id}_fastp_R2.fastq.gz \
        -j ${sample_id}_fastp.json \
        -h ${sample_id}_fastp.html \
        -q 30 \
        -e 25 \
        -n 5 \
        -l 40 \
        -c \
        -x \
        -p \
        --verbose \
        -w ${task.cpus}
    rm ${sample_id}_raw_R1.fastq.gz
    rm ${sample_id}_raw_R2.fastq.gz
    """  
}
