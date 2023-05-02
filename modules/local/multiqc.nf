
process multiqc {

    publishDir "${params.outDir}/multiqc", mode:'symlink'
   
    input:
        path("*")
    
    output:
        file("multiqc_all.html")

    script:
        """

        multiqc * \
            -c ${baseDir}/assets/multiqc_config.yaml \
            -n multiqc_all
            
	"""
}
