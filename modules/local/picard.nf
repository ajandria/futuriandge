
 process picard_metrics {

     publishDir "${params.outDir}/picard", mode:'symlink'

     tag "CollectRnaSeqMetrics on ${meta}"
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: picard_to_multiqc

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
             strandedness = meta.single_end == true ? 'SECOND_READ_TRANSCRIPTION_STRAND' : 'SECOND_READ_TRANSCRIPTION_STRAND'
         } else if (meta.strandedness == 'forward') {
             strandedness = meta.single_end == true ? 'FIRST_READ_TRANSCRIPTION_STRAND' : 'FIRST_READ_TRANSCRIPTION_STRAND'
         } else if (meta.strandedness == 'unstranded') {
             strandedness = meta.single_end == true ? 'NONE' : 'NONE'
         }
     } else {
         log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to 'NONE'."
     }

    """
    picard CollectRnaSeqMetrics -I ${bam} -O ${meta.id}_picard_CollectRnaSeqMetrics \
        --REF_FLAT ${params.refFlat} \
        --STRAND ${strandedness}
    """ 
 }
