
include { DESEQ2 } from "${baseDir}/modules/local/downstream_processing.nf"

workflow DOWNSTREAM_PROCESSING {

    take:
    reads

    main:
    DESEQ2(reads)

}
