
process multiqc {

    conda 'conda-forge::importlib-metadata bioconda::multiqc=1.14'

    publishDir "${params.outDir}/multiqc", mode:'symlink'
   
    input:
        path("*")
        path("*")
        path("*")
        path("*")
        path("*")
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
