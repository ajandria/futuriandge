 
 process featureCounts {

     publishDir "${params.outDir}/featureCounts", mode:'symlink'
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: featureCounts_to_multiqc

        script:
        """
         featureCounts -t gene \
            -s 2 \
            -T ${task.cpus} \
            --verbose \
            -p \
            -a ${params.gtf} \
            -o ${meta.id}_featureCounts_matrix.txt \
            ${bam}
	    """
 }
