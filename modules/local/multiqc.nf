
process multiqc {

    publishDir "${params.outDir}/multiqc", mode:'symlink'
   
    input:
        path("*")
    
    output:
        file("multiqc_all.html")

    script:
        """
        multiqc * \
            -n multiqc_all

	"""
}
