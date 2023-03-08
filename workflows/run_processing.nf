
include { run_pipeline as RUN_PROCESSING } from "${baseDir}/subworkflows/run_pipeline.nf"
include { run_MultiQC as MULTIQC } from "${baseDir}/subworkflows/run_MultiQC.nf"

workflow RUN_PROCESSING_FOR_DE_ANALYSIS {

    // Processing workflow
    RUN_PROCESSING(
        // 0. Read paired end files
        Channel
            .fromFilePairs(params.reads, checkIfExists: true)
    )

    // Gather QC stats
    MULTIQC(
        RUN_PROCESSING.out.qc_MultiQC_input
    )
}
