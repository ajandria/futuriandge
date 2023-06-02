#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 ajandria/futuriandge - Bulk RNA-seq pipeline for differential gene expression analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ajandria/futuriandge
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// Initial params
params.contrasts = null // if null dont run DESeq2 downstream subworkflow
params.input = "${baseDir}/tests/test_sample_sheet.csv"
params.outDir = "${baseDir}/results"
params.organism = null

// Constant for rRNA screening using SILVA dbs
params.rRNA_db_path = '/archive/users/ajan/references/rRNA_databases'

if (params.organism == null) {
    error "params.organism was not provided. Please specify an organism."
} else if (params.organism == 'human') {
    println("Organism is human")
    params.reference_genome = '/archive/users/ajan/references/Homo_sapiens_GRCh38_109/Homo_sapiens.GRCh38.109_genomeDir_STAR_2.7.10b'
    params.gtf = '/archive/users/ajan/references/Homo_sapiens_GRCh38_109/Homo_sapiens.GRCh38.109.chr.gtf'
    params.refFlat = '/archive/users/ajan/references/Homo_sapiens_GRCh38_109/Homo_sapiens.GRCh38.109.chr.gtf.gz.refflat.gz'
    params.gene_map = '/archive/users/ajan/references/Homo_sapiens_GRCh38_109/homo_sapiens_salmon_1.10.1/txp2gene.tsv'
    params.salmon_index = '/archive/users/ajan/references/Homo_sapiens_GRCh38_109/homo_sapiens_salmon_1.10.1/salmon_index-with_decoys'
} else if (params.organism == 'mouse') {
    println("Organism is mouse")
    params.reference_genome = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108_genomeDir_STAR_2.7.10b'
    params.gtf = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf'
    params.refFlat = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf.refFlat.gz'
    params.gene_map = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/salmon_1.10_1/txp2gene.tsv'
    params.salmon_index = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/salmon_1.10_1/salmon_index'
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)


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
| Sample sheet                         : $params.samplesheet
|
=====================================================================================================================
""".stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UNIFORMAL } from "${baseDir}/workflows/run_processing.nf"

//
// WORKFLOW: Run main ajandria/futuriandge analysis pipeline
//
workflow FUTURIANDGE {
    UNIFORMAL()
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
    FUTURIANDGE()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/