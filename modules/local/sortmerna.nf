// Define the rRNA database paths
rRNA = [
    "rfam-5.8s-database-id98.fasta",
    "rfam-5s-database-id98.fasta",
    "silva-arc-16s-id95.fasta",
    "silva-arc-23s-id98.fasta",
    "silva-bac-16s-id90.fasta",
    "silva-bac-23s-id98.fasta",
    "silva-euk-18s-id95.fasta",
    "silva-euk-28s-id98.fasta"
]

process sortmerna {
    tag "sortmerna on ${meta.id}"  // Tag the process with the id of the input data item
    label 'process_high'  // Resources are configured in the base config
    publishDir "${params.outDir}/sortmerna", mode:'symlink'  // Specify the output directory

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.4--h9ee0642_0' :
        'biocontainers/sortmerna:4.3.4--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)  // Accept a tuple of metadata and one or two FASTQ files

    output:
    tuple val(meta), path("${meta.id}_*.fq.gz"), emit: reads  // Emit the id and the cleaned FASTQ file(s)
    file("*.log") // Emit the log files

    script:
    def input_reads = meta.single_end == 'true' ? "--reads ${reads[0]}" : "--reads ${reads[0]} --reads ${reads[1]}"
    
    // Prepare rRNA databases for the --ref option
    def refs = rRNA.collect { "--ref ${it}" }.join(' ')  // Collect reference databases into a string

    """
    echo "Reference databases: ${refs}"

    if [ -d ${params.outDir}/sortmerna/${meta.id} ]; then
    rm -r ${params.outDir}/sortmerna/${meta.id}
    fi

    mkdir -p ${params.outDir}/sortmerna/${meta.id}

    sortmerna \
    ${refs} \
    ${input_reads} \
    --workdir ${params.outDir}/sortmerna/${meta.id} \
    -threads ${task.cpus} \
    --fastx \
    -v \
    -m 31744 \
    --out2 \
    --paired_out \
    --other ${meta.id} \
    --aligned contamined_${meta.id}
    """
}
