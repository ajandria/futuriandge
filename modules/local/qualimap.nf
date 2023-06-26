process qualimap {
    tag "${meta.id}"
    label 'intense'

    publishDir "${params.outDir}/qualimap", mode: 'symlink'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.2.2d--1' :
        'biocontainers/qualimap:2.2.2d--1' }"
   
    input:
    tuple val(meta), path(bam)
    
    output:
    tuple val(meta), path("${meta.id}"), emit: qualimap_to_multiqc

    script:
    // Dict of library types
    def strandedness_opts = ['reverse', 'forward', 'unstranded']
    
    // By default, try to auto-detect library type
    def strandedness = 'non-strand-specific'

    // Check if library type is valid and assign library type
    if (strandedness_opts.contains(meta.strandedness)) {
        if (meta.strandedness == 'reverse') {
            strandedness = 'strand-specific-reverse'
        } else if (meta.strandedness == 'forward') {
            strandedness = 'strand-specific-forward'
        } else if (meta.strandedness == 'unstranded') {
            strandedness = 'non-strand-specific'
        }
    } else {
        log.info "[qualimap rnaseq] Invalid library type specified '--libType=${lib_type}', defaulting to 'non-strand-specific'."
    }

    // Check if data is paired-end or single-end for -pe arg
    def is_paired_end = !meta.single_end ? '-pe' : ''

    """
    mkdir -p ${params.outDir}/qualimap/${meta.id}

    export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap/${meta.id}"

    qualimap rnaseq \
       -bam ${bam} \
       -gtf ${params.gtf} \
       -outdir ${meta.id} \
       -p ${strandedness} \
       -outformat PDF:HTML \
       --java-mem-size=16G \
       ${is_paired_end}
    """
}
