
include { fastqc } from "${baseDir}/modules/local/fastqc.nf"
include { fastp } from "${baseDir}/modules/local/fastp.nf"
include { star } from "${baseDir}/modules/local/star.nf"
include { qualimap } from "${baseDir}/modules/local/qualimap.nf"
include { picard_matrix } from "${baseDir}/modules/local/picard.nf"
include { samtools } from "${baseDir}/modules/local/samtools.nf"
include { featureCounts } from "${baseDir}/modules/local/subread.nf"

workflow run_pipeline {
    
    take:
    fastq_files

    main:
    // 1. FastQC on raw .fastq files
    fastqc(fastq_files)
    // 2. Run fastp on the raw .fastq files
    fastp(fastq_files)
    // 3. Align with STAR
    star(fastp.out[0])
    // 4. Stats with qualimap
    qualimap(star.out[0])
    // 5. Stats with picard
    picard_matrix(star.out[0])
    // 6. Stats with samtools
    samtools(star.out[0])
    // 7. Count with subread
    featureCounts(star.out[0])

    emit:
    qc_MultiQC_input = (featureCounts.out.featureCounts_to_multiqc.collect()
    .mix(samtools.out.samtools_to_multiqc.collect())
    .mix(picard_matrix.out.picard_to_multiqc.collect())
    .mix(qualimap.out.qualimap_to_multiqc.collect())
    .mix(star.out.star_to_multiqc.collect())
    .mix(fastp.out.fastp_to_multiqc.collect())
    .mix(fastqc.out.fastqc_on_raw_to_multiqc.collect())).collect()

}
