#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SodzieR/DE_RNA-Seq_nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/SodzieR/DE_RNA-Seq_nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.reads = "${baseDir}/FastQ/*_{R1,R2}.fastq.gz"
params.outDir = "${baseDir}/results"
params.organism = null
params.strandedness = null
params.protocol = null
params.reference_genome = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108_genomeDir_STAR_2.7.10b'
params.gtf = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf'
params.refFlat = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf.refFlat.gz'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info """\
====================================================================================================================

   U N I F O R M A L   B U L K   R N A - S E Q   D I F F E R E N T I A L   G E N E   
   E X P R E S S I O N   A N A L Y S I S   P I P E L I N E

   V e r s i o n:   0 . 0 . 1 . d


=====================================================================================================================
|
| .fastq files                          : $params.reads
|
=====================================================================================================================
""".stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RUN_PROCESSING_FOR_DE_ANALYSIS as RUN_PIPELINE } from "${baseDir}/workflows/run_processing.nf"

//
// WORKFLOW: Run main nf-core/rnaseq analysis pipeline
//
workflow DE_RNA_Seq_nf {
    RUN_PIPELINE()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    DE_RNA_Seq_nf()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/