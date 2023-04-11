// Modified based od nf-core/salmon_quant module
process SALMON_QUANT {

    publishDir "${params.outDir}/salmon_quant", mode:'symlink'
    
    tag "${sample_id}"
    label "intense"

    conda "conda-forge::boost-cpp bioconda::salmon=1.10.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salmon:1.10.1--h7e5ed60_0' :
        'quay.io/biocontainers/salmon:1.10.1--h7e5ed60_0' }"

    input:
    tuple val(sample_id), path(reads)

    output:
    path  ("${sample_id}_salmon")                          , emit: salmon_to_multiqc
    path  "${sample_id}_salmon_versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    if({params.strandedness == 'reverse'} && {params.protocol = 'paired-end'})

        """
            salmon quant \
                -i ${params.salmon_index} \
                -l ISR \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${params.gene_map}

            cat <<-END_VERSIONS > ${sample_id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS

	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'paired-end'})

        """
            salmon quant \
                -i ${salmon_index} \
                -l ISF \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${gene_map}

            cat <<-END_VERSIONS > ${sample_id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'paired-end'})

        """
            salmon quant \
                -i ${salmon_index} \
                -l IU \
                -1 ${reads[0]} \
                -2 ${reads[1]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${gene_map}

            cat <<-END_VERSIONS > ${sample_id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS

	    """

    else if({params.strandedness == 'reverse'} && {params.protocol = 'single-end'})

        """
            salmon quant \
                -i ${salmon_index} \
                -l SR \
                -1 ${reads[0]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${gene_map}

            cat <<-END_VERSIONS > ${sample_id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS
	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'single-end'})

        """
            salmon quant \
                -i ${salmon_index} \
                -l SF \
                -1 ${reads[0]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${gene_map}

            cat <<-END_VERSIONS > ${sample_id}_salmon_versions.yml
            "${task.process}":
                salmon: \$(echo \$(salmon --version) | sed -e "s/salmon //g")
            END_VERSIONS
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'single-end'})

        """
            salmon quant \
                -i ${salmon_index} \
                -l U \
                -1 ${reads[0]} \
                -p ${task.cpus} \
                -o ${sample_id}_salmon \
                -g ${gene_map}
	    """

    else

        throw new IllegalArgumentException("Unknown strandedness $params.strandedness")
        
 }
