#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SodzieR/DE_RNA-Seq_nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This pipeline has been based on the ENCODE bulk RNA-seq data processing standards
        https://www.encodeproject.org/data-standards/encode4-bulk-rna/
    Github : https://github.com/SodzieR/DE_RNA-Seq_nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.reads = "${baseDir}/FastQ/*_{R1,R2}.fastq.gz"
params.outDir = "${baseDir}/results"
params.organism = null
params.strandedness = null
params.protocol = null
params.reference_genome = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108_genomeDir_STAR_2.7.10b'
params.gtf = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf'
params.refFlat = '/archive/users/ajan/references/Mus_musculus.GRCm39.108/Mus_musculus.GRCm39.108.chr.gtf.refFlat.gz'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info """\
====================================================================================================================

   U N I F O R M A L   B U L K   R N A - S E Q   D I F F E R E N T I A L   G E N E   
   E X P R E S S I O N   A N A L Y S I S   P I P E L I N E

   V e r s i o n:   0 . 0 . 1 . d


=====================================================================================================================
|                                                                                                                   
| .fastq files                          : $params.reads                                                                                         
|                                                                                                                     
=====================================================================================================================
""".stripIndent()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { fastqc } from "${baseDir}/modules/local/fastqc.nf"
include { fastp } from "${baseDir}/modules/local/fastp.nf"
include { star } from "${baseDir}/modules/local/star.nf"
include { qualimap } from "${baseDir}/modules/local/qualimap.nf"
include { picard_matrix } from "${baseDir}/modules/local/picard.nf"
include { samtools } from "${baseDir}/modules/local/samtools.nf"
include { featureCounts } from "${baseDir}/modules/local/subread.nf"
include { multiqc } from "${baseDir}/modules/local/multiqc.nf"

workflow {

    // 0. Read paired end files
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set{read_pairs_ch}

    // 1. FastQC on raw .fastq files
    fastqc(read_pairs_ch)

    // 2. Run fastp on the raw .fastq files
    fastp(read_pairs_ch)

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

    // 8. Run MultiQC an aggregate all of the stats for samples
    multiqc_ch = multiqc(fastqc.out[1].collect(),
        fastp.out[1].collect(),
        star.out[1].collect(),
        qualimap.out[0].collect(),
        picard_matrix.out[0].collect(),
        samtools.out[0].collect(),
        featureCounts.out[0].collect()
        )

}

