
include { DESeq2_DOWNSTREAM } from "${baseDir}/modules/local/downstream_processing.nf"

workflow DOWNSTREAM_PROCESSING {

    take:
    reads

    main:
    DESeq2_DOWNSTREAM(reads)

}
