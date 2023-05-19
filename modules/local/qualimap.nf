
 process qualimap {

    tag "${meta}"

    label = 'intense'

     publishDir "${params.outDir}/qualimap", mode:'symlink'
   
     input:
         tuple val(meta), path(bam)
    
     output:
         path("*"), emit: qualimap_to_multiqc

    script:
     // Dict of library types
     def strandedness_opts = [
         'reverse', 'forward', 'unstranded'
     ]
    
     // By default try to auto-detect library type
     def lib_type =  'non-strand-specific'

     // Check if library type is valid and assign library type
     if ({strandedness_opts.contains(meta.strandedness)}) {
         if (meta.strandedness == 'reverse') {
             strandedness = meta.single_end == true ? 'strand-specific-reverse' : 'strand-specific-reverse'
         } else if (meta.strandedness == 'forward') {
             strandedness = meta.single_end == true ? 'strand-specific-forward' : 'strand-specific-forward'
         } else if (meta.strandedness == 'unstranded') {
             strandedness = meta.single_end == true ? 'non-strand-specific' : 'non-strand-specific'
         }
     } else {
         log.info "[Salmon Quant] Invalid library type specified '--libType=${lib_type}', defaulting to 'non-strand-specific'."
     }
        """
        mkdir -p ${params.outDir}/qualimap
         export JAVA_OPTS="-Djava.io.tmpdir=${params.outDir}/qualimap"
         qualimap rnaseq \
            -bam ${bam} \
            -gtf ${params.gtf} \
            -outdir ${meta.id} \
            -p ${strandedness} \
            -pe \
            -outformat PDF:HTML \
            --java-mem-size=16G
	    """        
 }
