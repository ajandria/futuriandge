
include { INPUT_CHECK } from "${baseDir}/subworkflows/local/input_check.nf"
include { PROCESSING } from "${baseDir}/subworkflows/local/run_pipeline.nf"
include { DOWNSTREAM_PROCESSING } from "${baseDir}/subworkflows/local/downstream_analysis.nf"
include { QC_MULTIQC } from "${baseDir}/subworkflows/local/run_MultiQC.nf"

workflow UNIFORMAL {

    // Check samplesheet
    INPUT_CHECK(
        params.input
    )

    // Processing workflow
    PROCESSING(
        // 0. Read paired end files
        INPUT_CHECK.out.reads
    )

    // Gather QC stats
    QC_MULTIQC(
        PROCESSING.out.qc_MultiQC_input
    )

    // If metadata is provided, proceed to downstream analysis
    if ( params.contrasts != null ) {
        DOWNSTREAM_PROCESSING(
            PROCESSING.out.featureCounts_paths
        )
    }

}
