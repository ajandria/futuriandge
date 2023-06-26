process FEATURECOUNTS {
    publishDir "${params.outDir}/featureCounts", mode: 'symlink'
    
    tag "featureCounts on ${meta.id}" 

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/subread:2.0.1--hed695b0_0' :
        'biocontainers/subread:2.0.1--hed695b0_0' }"

    input:
    tuple val(meta), path(bam)
    
    output:
    path("*"), emit: featureCounts_to_multiqc
    path("${meta.id}_downstream_info.txt"), emit: counts_path
    
    script:
    def strandedness_opts = ['reverse', 'forward', 'unstranded']
    def strandedness = '0'  // By default try to auto-detect library type
    
    if (strandedness_opts.contains(meta.strandedness)) {
        if (meta.strandedness == 'reverse') {
            strandedness = meta.single_end ? '2' : '2'  // Reverse stranded
        } else if (meta.strandedness == 'forward') {
            strandedness = meta.single_end ? '1' : '1'  // Forward stranded
        } else if (meta.strandedness == 'unstranded') {
            strandedness = meta.single_end ? '0' : '0'  // Unstranded
        }
    } else {
        log.info "[subread featureCounts] Invalid library type specified '--libType=${strandedness}', defaulting to '0'."
    }
    
    """
    echo "Strandedness: ${strandedness}"

    featureCounts -t gene \
        -s ${strandedness} \
        -T ${task.cpus} \
        --verbose \
        -p \
        -a ${params.gtf} \
        -o ${meta.id}_featureCounts_matrix.txt \
        ${bam}
        
    # Define column names
    echo "sample_id,is_single_end,strandedness,count_matrix_path" > ${meta.id}_downstream_info.txt
    
    # Write sample info to file
    echo "${meta.id},${meta.single_end},${meta.strandedness},${PWD}/${meta.id}_featureCounts_matrix.txt" >> ${meta.id}_downstream_info.txt
    """
}
