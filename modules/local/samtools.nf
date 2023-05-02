 
  process samtools {

    tag "${meta}"

     publishDir "${params.outDir}/samtools", mode:'symlink'
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: samtools_to_multiqc

     script:
         """
         samtools flagstat ${bam} -@ 4 > ${meta.id}_flagstat.txt
	"""
 }
