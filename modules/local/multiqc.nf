
process MULTIQC {

    publishDir "${params.outDir}/multiqc", mode:'symlink'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.14--pyhdfd78af_0' :
        'biocontainers/multiqc:1.14--pyhdfd78af_0' }"
   
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
