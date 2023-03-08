 
  process samtools {

    conda 'conda-forge::zlib bioconda::samtools=1.16.1'

     publishDir "${params.outDir}/samtools", mode:'symlink'
   
     input:
         tuple val(sample_id), path(bam)
    
     output:
         path("*"), emit: samtools_to_multiqc

     script:
         """
         samtools flagstat ${bam} -@ 4 > ${sample_id}_flagstat.txt
	"""
 }
