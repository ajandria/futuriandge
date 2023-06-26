
process MULTIQC {

    publishDir "${params.outDir}/multiqc", mode:'symlink'
   
    input:
        path("*")
    
    output:
        file("multiqc_*.html")

    script:
        """
        multiqc * \
            -c ${baseDir}/assets/multiqc_config.yaml \
            -n multiqc_all

        multiqc * \
            -c ${baseDir}/assets/multiqc_config.yaml \
            -n multiqc_all_interactive \
            --interactive
	"""
}
