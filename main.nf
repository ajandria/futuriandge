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
params.reads = "${baseDir}/FastQ/*_{R1,R2}.fastq*"
params.outDir = "${baseDir}/results"
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

process fastqc {

    publishDir "${params.outDir}/FastQC", mode: 'symlink'

    conda 'conda-forge::openjdk bioconda::fastqc'

    tag "FastQC on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*"), emit: fastqcRaw_to_multiqc

    """
    fastqc -t ${task.cpus} ${reads}
    """  
}

process fastp {

    publishDir "${params.outDir}/fastp", mode: 'symlink'

    conda 'bioconda::fastp'

    tag "fastp on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_R*"), emit: fastp_to_align
    path("*"), emit: fastp_to_multiqc

    """
    # Run fastp
    fastp -i ${reads[0]} -I ${reads[1]} -o ${sample_id}_R1.fastq.gz -O ${sample_id}_R2.fastq.gz \
        -j ${sample_id}_fastp.json \
        -h ${sample_id}_fastp.html \
        -q 30 \
        -e 25 \
        -n 5 \
        -l 40 \
        -c \
        -x \
        -w ${task.cpus}

    # Rename R1 to just sample_name
    #sed -i '' 's/_R1//g' ${sample_id}_fastp.json # for MacOS (due to backup file creation)
    #sed -i 's/_R1//g' ${sample_id}_fastp.json # for Linux
    """  
}

process star {

    publishDir "${params.outDir}/star", mode: 'symlink'

    conda = 'bioconda::star=2.1.10b'

    tag "star on ${sample_id}"

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*")

    """
    STAR --runMode alignReads \
    --genomeDir $star_index_in \
    --sjdbGTFfile $gtf_in \
    --readFilesIn ${reads[0]} ${reads[1]} \
    # IF GZIPPED # --readFilesCommand zcat \
    --readStrand Reverse \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ${sample_id} \
    --runThreadN ${task.cpus}
    """  
}

 process qualimap {

    conda 'bioconda::qualimap'

     publishDir "${params.outDir}/qualimap", mode:'copy'
   
     input:
         tuple val(sample_id), file(bam) from star_qualimap
    
     output:
         path id into qualimap_out_path

     script:
         """
         export JAVA_OPTS="-Djava.io.tmpdir=/home/ajan/shared_folder/nextflow_liam_star/trash"
         qualimap rnaseq -bam ${bam} -gtf ${gtf_in} -outdir ${id}/ -pe -p strand-specific-reverse -outformat PDF:HTML --java-mem-size=16G
	"""
 }

 process picard_matrix {

    conda 'bioconda::picard'

     publishDir "${params.outDir}/picard", mode:'copy'
   
     input:
         set val(id), file(bam) from star_picard
    
     output:
         path id into picard_out_path

     script:
         """
         picard CollectRnaSeqMetrics -I ${bam} -O ${id} --REF_FLAT ${ref_flat_in} --STRAND SECOND_READ_TRANSCRIPTION_STRAND
	"""
 }

 process samtools_index {

    conda 'bioconda::samtools'

     publishDir "${params.outDir}/star", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_index
    
     output:
         file("*")

     script:
         """
         samtools index -@ 4 ${bam}
	"""
 }

 process samtools_flagstat {

    conda 'bioconda::samtools'

     publishDir "${params.outDir}/samtools", mode:'copy'
   
     input:
         set val(id), file(bam) from star_samtools_flagstat
    
     output:
         path ("${id}.txt") into samtools_flagstat_out_path

     script:
         """
         samtools flagstat ${bam} -@ 4 > ${id}.txt
	"""
 }

 process featureCounts {

    conda 'bioconda::subread'

     publishDir "${params.outDir}/featureCounts", mode:'copy'
   
     input:
         set val(id), file(bam) from star_featureCounts
    
     output:
         path ("*.summary") into featureCounts_multiqc
         file ("*")

     script:
         """
         featureCounts -t gene \
            -s 2 \
            -T ${task.cpus} \
            --verbose \
            -p \
            -a ${gtf_in} \
            -o ${id}_featureCounts_matrix.txt \
            ${bam}
	"""
 }


process multiqc {

    conda 'conda-forge::importlib-metadata bioconda::multiqc'

    publishDir "${params.outDir}/multiqc", mode:'symlink'
   
    input:
        path("*")
        path("*")
    
    output:
        file("multiqc_all.html")

    script:
        """
        multiqc * \
            -n multiqc_all

	"""
}


workflow {

    // 0. Read paired end files
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set{read_pairs_ch}

    // 1. FastQC on raw .fastq files
    fastqc_ch = fastqc(read_pairs_ch)
    //fastqc_ch.view()

    // 2. Run fastp on the raw .fastq files
    fastp_ch = fastp(read_pairs_ch)
    //fastp_ch.inspect()

    // 3. Align with STAR
    star_ch = star(fastp_ch.fastp_to_align)
    star_ch.view()

    // X. Run MultiQC an aggregate all of the stats for samples
    //fastp_ch.fastp_to_multiqc.view()
    multiqc_ch = multiqc(fastqc_ch.fastqcRaw_to_multiqc.collect(),
        fastp_ch.fastp_to_multiqc.collect())

}

