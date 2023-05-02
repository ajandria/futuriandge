
include { INPUT_CHECK } from "${baseDir}/subworkflows/local/input_check.nf"
include { PROCESSING } from "${baseDir}/subworkflows/run_pipeline.nf"
include { MULTIQC } from "${baseDir}/subworkflows/run_MultiQC.nf"

workflow UNIFORMAL {

    // Check samplesheet
    INPUT_CHECK(
        params.samplesheet
    )

    // Processing workflow
    PROCESSING(
        // 0. Read paired end files
        INPUT_CHECK.out.reads
    )

    // Gather QC stats
    MULTIQC(
        PROCESSING.out.qc_MultiQC_input
    )
}
