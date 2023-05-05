
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

    // Define input reads
    def input_reads = meta.single_end == 'true' ? "-i ${reads[0]}" : "-i ${reads[0]} -I ${reads[1]}"
    def output_reads = meta.single_end == 'true' ? "-o ${meta.id}_fastp_R1.fastq.gz" : "-o ${meta.id}_fastp_R1.fastq.gz -O ${meta.id}_fastp_R2.fastq.gz"

    // Dict of library types
    def strandedness_opts = [
        'reverse', 'forward', 'unstranded'
    ]
    
    // By default try to auto-detect library type
    def lib_type =  'A'

    // Check if library type is valid and assign library type
    if ({strandedness_opts.contains(meta.strandedness)}) {
        if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end == true ? 'ISR' : 'SR'
        } else if (meta.strandedness == 'forward') {
            strandedness = meta.single_end == true ? 'ISF' : 'SF'
        } else if (meta.strandedness == 'unstranded') {
            strandedness = meta.single_end == true ? 'IU' : 'U'
        }
    } else {
        log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
    }

        """
        cp ${reads[0]} ${meta.id}_raw_R1.fastq.gz
        cp ${reads[1]} ${meta.id}_raw_R2.fastq.gz

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
        rm ${meta.id}_raw_R2.fastq.gz
        """
}