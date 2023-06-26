
include { MULTIQC } from "${baseDir}/modules/local/multiqc.nf"

workflow QC_MULTIQC {

    take:
    input_for_MultiQC

    main:
    MULTIQC(input_for_MultiQC)

}
