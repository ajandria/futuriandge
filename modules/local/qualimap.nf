
 process qualimap {

    label = 'intense'

     publishDir "${params.outDir}/qualimap", mode:'symlink'
   
     input:
         tuple val(sample_id), path(bam)
    
     output:
         path("*"), emit: qualimap_to_multiqc

    script:
    if({params.strandedness == 'reverse'} && {params.protocol = 'paired-end'})

        """
        mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p strand-specific-reverse \
            -pe \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'paired-end'})

        """
         mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p strand-specific-forward \
            -pe \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'paired-end'})

        """
         mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p non-strand-specific \
            -pe \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else if({params.strandedness == 'reverse'} && {params.protocol = 'single-end'})

        """
         mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p strand-specific-reverse \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else if({params.strandedness == 'forward'} && {params.protocol = 'single-end'})

        """
         mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p strand-specific-forward \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else if({params.strandedness == 'unstranded'} && {params.protocol = 'single-end'})

        """
         mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${sample_id} \
            -p non-strand-specific \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """

    else

        throw new IllegalArgumentException("Unknown strandedness $params.strandedness")
        
 }
