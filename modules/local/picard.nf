
 process picard_matrix {

     publishDir "${params.outDir}/picard", mode:'symlink'

     tag "qualimap on ${sample_id}"
   
     input:
         tuple val(sample_id), path(bam)
    
     output:
         path("*"), emit: picard_to_multiqc

     script:
         """
         picard CollectRnaSeqMetrics -I ${bam} -O ${sample_id}_picard_CollectRnaSeqMetrics --REF_FLAT ${params.refFlat} --STRAND SECOND_READ_TRANSCRIPTION_STRAND
	"""
 }
