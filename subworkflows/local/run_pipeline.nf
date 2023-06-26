
include { FASTQC } from "${baseDir}/modules/nf-core/fastqc/main.nf"
include { FASTP } from "${baseDir}/modules/local/fastp.nf"
include { SORTMERNA } from "${baseDir}/modules/local/sortmerna.nf"
include { STAR_ALIGN } from "${baseDir}/modules/local/star.nf"
include { QUALIMAP } from "${baseDir}/modules/local/qualimap.nf"
include { PICARD_METRICS } from "${baseDir}/modules/local/picard.nf"
include { SAMTOOLS_FLAGSTAT } from "${baseDir}/modules/local/samtools.nf"
include { FEATURECOUNTS } from "${baseDir}/modules/local/subread.nf"
include { SALMON_QUANT } from "${baseDir}/modules/nf-core/salmon/quant/main.nf"

workflow PROCESSING {
    
    take:
    fastq_files

    main:
    // 1. FastQC on raw .fastq files
    FASTQC(fastq_files)
    // 2. Run fastp on the raw .fastq files
    FASTP(fastq_files)
    // 2. Remove rRNA
    SORTMERNA(FASTP.out[0])
    // 3. Align with STAR
    STAR_ALIGN(SORTMERNA.out[0])
    // 4. Stats with qualimap
    QUALIMAP(STAR_ALIGN.out[0])
    // 5. Stats with picard
    PICARD_METRICS(STAR_ALIGN.out[0])
    // 6. Stats with samtools
    SAMTOOLS_FLAGSTAT(STAR_ALIGN.out[0])
    // 7. Count with subread
    FEATURECOUNTS(STAR_ALIGN.out[0])
    // 8. Additionally quantify with salmon
    SALMON_QUANT(SORTMERNA.out[0])

    emit:
    qc_MultiQC_input = (FEATURECOUNTS.out.featureCounts_to_multiqc.collect()
    .mix(SAMTOOLS_FLAGSTAT.out.samtools_to_multiqc.map { it -> it[1] }.collect())
    .mix(PICARD_METRICS.out.picard_to_multiqc.map { it -> it[1] }.collect())
    .mix(QUALIMAP.out.qualimap_to_multiqc.map { it -> it[1] }.collect())
    .mix(STAR_ALIGN.out.star_to_multiqc.collect())
    .mix(SORTMERNA.out[1].collect())
    .mix(FASTP.out.fastp_json.collect())
    .mix(FASTQC.out.zip.map { it -> it[1] }.collect())
    .mix(SALMON_QUANT.out.salmon_to_multiqc.collect())).collect()

    //featureCounts_meta = featureCounts.out.counts_path.map { it -> it[0] }.collect()
    //featureCounts_paths = featureCounts.out.counts_path.map { it -> it[1] }.collect()
    featureCounts_paths = FEATURECOUNTS.out.counts_path.collect()

}
