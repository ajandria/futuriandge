
 process picard_metrics {

     publishDir "${params.outDir}/picard", mode:'symlink'

     tag "CollectRnaSeqMetrics on ${meta}"
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: picard_to_multiqc

     script:
        """
        picard CollectRnaSeqMetrics -I ${bam} -O ${meta.id}_picard_CollectRnaSeqMetrics \
            --REF_FLAT ${params.refFlat} \
            --STRAND SECOND_READ_TRANSCRIPTION_STRAND
	    """ 
 }
