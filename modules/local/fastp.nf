process fastp {

    // Tagging the process for easier identification in logs
    tag "fastp on ${meta.id}"

    // Define the intensity of resources required by the process
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    // Define the directory to publish results
    publishDir "${params.outDir}/fastp", mode: 'symlink'

    input:
    // Receive a tuple containing metadata and the path of input reads
    tuple val(meta), path(reads)

    output:
    // Emit two types of output - one for alignment and one for multiqc
    tuple val(meta), path("${meta.id}_fastp_R*.fastq.gz"), emit: fastp_to_align
    path("*_fastp.json"), emit: fastp_json
    path("*_fastp.html"), emit: fastp_html

    script:
    // Define input and output reads
    // Single-end or paired-end reads are accommodated
    def input_reads = meta.single_end == 'true' ? "-i ${reads[0]}" : "-i ${reads[0]} -I ${reads[1]}"
    def output_reads = meta.single_end == 'true' ? "-o ${meta.id}_fastp_R1.fastq.gz" : "-o ${meta.id}_fastp_R1.fastq.gz -O ${meta.id}_fastp_R2.fastq.gz"

    // Main script
    // Copy input reads, run fastp with specified parameters, and remove the copied input reads
    """
    cp ${reads[0]} ${meta.id}_raw_R1.fastq.gz
    if [[ "${meta.single_end}" != 'true' ]]; then
        cp ${reads[1]} ${meta.id}_raw_R2.fastq.gz
    fi

    fastp ${input_reads} \
        ${output_reads} \
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
    if [[ "${meta.single_end}" != 'true' ]]; then
        rm ${meta.id}_raw_R2.fastq.gz
    fi
    """
}
