// Modified based od nf-core/salmon_quant module
process SALMON_QUANT {

    publishDir "${params.outDir}/salmon_quant", mode:'symlink'
    
    tag "${meta}"
    label "intense"

    input:
    tuple val(meta), path(reads)

    output:
    path  ("${meta.id}_salmon")                          , emit: salmon_to_multiqc
    path  "${meta.id}_salmon_versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // Define input reads
    def input_reads = meta.single_end ? "-r ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"

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
