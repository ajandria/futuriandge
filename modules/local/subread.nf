 
 process featureCounts {

     publishDir "${params.outDir}/featureCounts", mode:'symlink'
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: featureCounts_to_multiqc
         path("${meta.id}_downstream_info.txt"), emit: counts_path

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

         # Define column names
         echo "sample_id,is_single_end,strandedness,count_matrix_path" > ${meta.id}_downstream_info.txt

         # Write sample info to file
         echo "${meta.id},${meta.single_end},${meta.strandedness},${PWD}/${meta.id}_featureCounts_matrix.txt" >> ${meta.id}_downstream_info.txt
	    """
 }
