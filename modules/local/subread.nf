 
 process featureCounts {

     publishDir "${params.outDir}/featureCounts", mode:'symlink'
   
     input:
         tuple val(sample_id), path(bam)
    
     output:
         path("*"), emit: featureCounts_to_multiqc

        script:
    if({params.strandedness == 'reverse'} && {params.protocol = 'paired-end'})

        """
         featureCounts -t gene \
            -s 2 \
            -T ${task.cpus} \
            --verbose \
            -p \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'paired-end'})

        """
         featureCounts -t gene \
            -s 1 \
            -T ${task.cpus} \
            --verbose \
            -p \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'paired-end'})

        """
         featureCounts -t gene \
            -s 0 \
            -T ${task.cpus} \
            --verbose \
            -p \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else if({params.strandedness == 'reverse'} && {params.protocol = 'single-end'})

        """
         featureCounts -t gene \
            -s 2 \
            -T ${task.cpus} \
            --verbose \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'single-end'})

        """
         featureCounts -t gene \
            -s 1 \
            -T ${task.cpus} \
            --verbose \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'single-end'})

        """
         featureCounts -t gene \
            -s 0 \
            -T ${task.cpus} \
            --verbose \
            -a ${params.gtf} \
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	    """

    else

        throw new IllegalArgumentException("Unknown strandedness $params.strandedness")
        
 }
