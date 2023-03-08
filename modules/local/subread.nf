 
 process featureCounts {

    conda 'conda-forge::zlib bioconda::subread=2.0.3'

     publishDir "${params.outDir}/featureCounts", mode:'symlink'
   
     input:
         tuple val(sample_id), path(bam)
    
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
            -o ${sample_id}_featureCounts_matrix.txt \
            ${bam}
	"""
 }
 