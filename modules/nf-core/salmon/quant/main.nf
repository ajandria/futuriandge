// Modified based od nf-core/salmon_quant module
process SALMON_QUANT {

    publishDir "${params.outDir}/salmon_quant", mode:'symlink'
    
    tag "${meta}"
    label "intense"

    conda "conda-forge::boost-cpp bioconda::salmon=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    path  ("${meta.id}_salmon")                          , emit: salmon_to_multiqc
    path  "${meta.id}_salmon_versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Define input reads
    def input_reads = params.protocol == 'single-end' ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

    // Dict of library types
    def strandedness_opts = [
        'paired-end', 'single-end', 'reverse', 'forward', 'unstranded'
    ]
    
    // By default try to auto-detect library type
    def strandedness =  'A'

    // Check if library type is valid and assign library type
    if ({params.strandedness} && {params.protocol}) {
        if ({strandedness_opts.contains(params.strandedness)} && {strandedness_opts.contains(params.protocol)}) {
            if (params.strandedness == 'reverse') {
                strandedness = params.protocol == 'paired-end' ? 'ISR' : 'SR'
            } else if (params.strandedness == 'forward') {
                strandedness = params.protocol == 'paired-end' ? 'ISF' : 'SF'
            } else if (params.strandedness == 'unstranded') {
                strandedness = params.protocol == 'paired-end' ? 'IU' : 'U'
            }
        } else {
            log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to auto-detection with '--libType=A'."
        }
    }
        """
            salmon quant \
                -i ${params.salmon_index} \
                -l ${strandedness} \
                ${input_reads} \
                -p ${task.cpus} \
                -o ${meta.id}_salmon \
                -g ${params.gene_map}

            cat <<-END_VERSIONS > ${meta.id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS

	    """
 }
