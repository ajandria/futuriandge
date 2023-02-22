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
params.reads = "${baseDir}/../FastQ/*_{R1,R2}.fastq*"
params.outDir = "${baseDir}/../results"
params.organism = null
params.strandedness = null
params.protocl = null

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    EXECUTE WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

log.info """\
====================================================================================================================

   U N I F O R M A L   B U L K   R N A - S E Q   D I F F E R E N T I A L   G E N E   
   E X P R E S S I O N   A N A L Y S I S   P I P E L I N E

   V e r s i o n:   0 . 0 . 1


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

process fastqc {

    publishDir "${params.outDir}/FastQC", mode: 'symlink'

    conda 'conda-forge::openjdk bioconda::fastqc'

    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*")

    script:
    """
    fastqc -t ${task.cpus} $reads
    """  
}

process fastp {

    publishDir "${params.outDir}/fastp", mode: 'symlink'

    conda 'bioconda::fastp'

    tag "fastp on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R*")
    path("*")

    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_R1_fastq.gz -O ${sample_id}_R2.fastq.gz -j ${sample_id}.json -h ${sample_id}.html \
        -q 30 \
        -e 25 \
        -n 5 \
        -l 40 \
        -c \
        -x \
        -w ${task.cpus}
    """  
}

process star {

    publishDir "${params.outDir}/star", mode: 'symlink'

    conda = '/home/ajan/.conda/envs/STAR_2.7.0d'

    tag "star on ${sample_id}"

    input:
    tuple val(sample_id), files(reads)

    output:
    path("*")

    script:
    """
    STAR --runThreadN ${task.cpus} --runMode alignReads \
    --genomeDir $star_index_in \
    --readFilesIn ${reads[0]} ${reads[1]} \
    --readStrand Reverse \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile $gtf_in \
    --sjdbOverhang 100 \
    --outFileNamePrefix ${id}
    """  
}

workflow {

    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set{read_pairs_ch}

    fastqc_ch = fastqc(read_pairs_ch)
    fastp(read_pairs_ch))
    //star_ch = star(fastp_ch)
}

