
 process qualimap {

    label = 'intense'

     publishDir "${params.outDir}/qualimap", mode:'symlink'
   
     input:
         tuple val(sample_id), path(bam)
    
     output:
         path("*"), emit: qualimap_to_multiqc

     script:
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
 }
