
include { MULTIQC } from "${baseDir}/modules/local/multiqc.nf"

workflow MULTIQC {

    take:
    input_for_MultiQC

    main:
    multiqc(input_for_MultiQC)

}
