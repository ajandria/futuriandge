
process fastqc {

    publishDir "${params.outDir}/FastQC", mode: 'symlink'

    conda 'conda-forge::openjdk bioconda::fastqc=0.12.1'

    label = 'intense'

    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*")

    """
    fastqc -t ${task.cpus} ${reads}
    """  
}
