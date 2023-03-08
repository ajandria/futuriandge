
include { multiqc } from "${baseDir}/modules/local/multiqc.nf"

workflow run_MultiQC {

    take:
    input_for_MultiQC

    main:
    multiqc(input_for_MultiQC)

}
