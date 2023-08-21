
process DESEQ2 {

    publishDir "${params.outDir}/DESeq2_DOWNSTREAM", mode: 'symlink'

    label = 'intense'

    input:
    path("*")

    output:
    path("*"), emit: DESeq2_DOWNSTREAM_out

    script:
    def org = params.organism
    def user_file = params.contrasts
    """
    # Merge all downstream info files
    cat *_downstream_info.txt > all_downstream_info.txt

    # Remove lines containing the header except for the first occurrence in the first line
    sed -i '1!{/sample_id,is_single_end,strandedness,count_matrix_path/d}' all_downstream_info.txt

    index-downstream.R all_downstream_info.txt ${org} ${user_file}

    """  
}
