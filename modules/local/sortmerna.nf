 
rRNA_1 = "${params.rRNA_db_path}/rfam-5.8s-database-id98.fasta"
rRNA_2 = "${params.rRNA_db_path}/rfam-5s-database-id98.fasta"
rRNA_3 = "${params.rRNA_db_path}/silva-arc-16s-id95.fasta"
rRNA_4 = "${params.rRNA_db_path}/silva-arc-23s-id98.fasta"
rRNA_5 = "${params.rRNA_db_path}/silva-bac-16s-id90.fasta"
rRNA_6 = "${params.rRNA_db_path}/silva-bac-23s-id98.fasta"
rRNA_7 = "${params.rRNA_db_path}/silva-euk-18s-id95.fasta"
rRNA_8 = "${params.rRNA_db_path}/silva-euk-28s-id98.fasta"

 process sortmerna {

     publishDir "${params.outDir}/sortmerna", mode:'symlink'
     
     label 'intense'

     tag "sortmerna on ${sample_id}"
 
     input:
         tuple val(sample_id), path(reads)
  
     output:
         tuple val(sample_id), path("${sample_id}_*.fq.gz")
         file("*.log")

     // --paired_out -> says to keep pairs together (if one matches, then count both as matching db) - no need for repair.sh
     // fwd and rev files are just the naming scheme used (in reality these are the R1 and R2 files)

    script:
    if(params.protocol == 'paired-end')

        """
         mkdir -p ${params.outDir}/sortmerna/${sample_id}
         sortmerna \
         --ref ${rRNA_1} \
         --ref ${rRNA_2} \
         --ref ${rRNA_3} \
         --ref ${rRNA_4} \
         --ref ${rRNA_5} \
         --ref ${rRNA_6} \
         --ref ${rRNA_7} \
         --ref ${rRNA_8} \
         --reads ${reads[0]} \
         --reads ${reads[1]} \
         --workdir ${params.outDir}/sortmerna/${sample_id} \
         -a ${task.cpus} \
         --fastx \
         -v \
         -m 31744 \
         --out2 \
         --paired_out \
         --other ${sample_id} \
         --aligned contamined_${sample_id}
	    """

    else if(params.protocol == 'single-end')

        """
         mkdir -p ${params.outDir}/sortmerna/${sample_id}
         sortmerna \
         --ref ${rRNA_1} \
         --ref ${rRNA_2} \
         --ref ${rRNA_3} \
         --ref ${rRNA_4} \
         --ref ${rRNA_5} \
         --ref ${rRNA_6} \
         --ref ${rRNA_7} \
         --ref ${rRNA_8} \
         --reads ${reads[0]} \
         --reads ${reads[1]} \
         --workdir ${params.outDir}/sortmerna/${sample_id} \
         -a ${task.cpus} \
         --fastx \
         -v \
         -m 31744 \
         --out2 \
         --paired_out \
         --other ${sample_id} \
         --aligned contamined_${sample_id}
	    """

    else

        throw new IllegalArgumentException("Unknown strandedness $params.strandedness")
     
 }
