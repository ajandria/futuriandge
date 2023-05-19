 
 process featureCounts {

     publishDir "${params.outDir}/featureCounts", mode:'symlink'
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: featureCounts_to_multiqc
         path("${meta.id}_downstream_info.txt"), emit: counts_path

     script:
     // Dict of library types
     def strandedness_opts = [
         'reverse', 'forward', 'unstranded'
     ]
    
     // By default try to auto-detect library type
     def lib_type =  'NONE'

     // Check if library type is valid and assign library type
     if ({strandedness_opts.contains(meta.strandedness)}) {
         if (meta.strandedness == 'reverse') {
             strandedness = meta.single_end == true ? '2' : '2'
         } else if (meta.strandedness == 'forward') {
             strandedness = meta.single_end == true ? '1' : '1'
         } else if (meta.strandedness == 'unstranded') {
             strandedness = meta.single_end == true ? '0' : '0'
         }
     } else {
         log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to 'NONE'."
     }
        """
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
