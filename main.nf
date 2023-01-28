#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SodzieR/DE_RNA-Seq_nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This pipeline has been based on the ENCODE bulk RNA-seq data processing standards
        https://www.encodeproject.org/data-standards/encode4-bulk-rna/
    Github : https://github.com/SodzieR/DE_RNA-Seq_nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.reads = "${baseDir}/*_{R1,R2}.fastq.gz"
params.outDir = "results"
params.organism = null
params.strandedness = null
params.protocl = null

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info """\
====================================================================================================================

   C E N T E R   F O R   R E G E N E R A T I V E   A N D   I M M U N E R E G U L A T I O N   M U B
   U N I F O R M A L   B U L K   R N A - S E Q   P I P E L I N E | v.0.1

=====================================================================================================================
|                                                                                                                   
| .fastq files                          : $fastq_in
| Genome                                :                                                                                            
|                                                                                                                     
=====================================================================================================================
""".stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process FastQC {

    publishDir "${params.outDir}/FastQC", mode:'copy'

    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("fastqc_${sample_id}_logs")

    script:
    """
    fastqc -t ${task.cpus} ${fastq_1} ${fastq_2}
    """  
}



