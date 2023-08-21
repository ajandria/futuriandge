#!/usr/bin/env Rscript

# Setup -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

downstream_ch_file <- args[1]
organism_in <- args[2]
user_file <- args[3]

if (organism_in == 'mouse') {
    organism <- 'mus_musculus'
} else if (organism_in == 'human') {
    organism <- 'homo_sapiens'
} else {
    stop('Provided organism not supported...')
}

#downstream_ch_file <- "/futuriandgeDownstream/used-locally-for-testing/counts/extended_source-pipeline-out_meta_file-dockerfile.txt"
#organism <- 'mus_musculus'
#user_file <- "/futuriandgeDownstream/used-locally-for-testing/template-user-request.xlsx"

sample_metadata <- readxl::read_excel(user_file, sheet = 'metadata')
contrasts <- readxl::read_excel(user_file, sheet = 'comparisons')

# Load --------------------------------------------------------------------
library(futuriandgeDownstream)

# F1
print('Running futuriandgeDownstream::return_count_matrix')
counts <- futuriandgeDownstream::return_count_matrix(downstream_ch_file)

# F2
print('Running futuriandgeDownstream::remove_mirna')
counts_no_mirnas <- futuriandgeDownstream::remove_mirna(counts, organism)

# F3
print('Running futuriandgeDownstream::run_dge')
dge_results <- futuriandgeDownstream::run_dge(count_matrix_raw = counts_no_mirnas,
                       metadata_raw = sample_metadata,
                       comparisons = contrasts)

# save dge results
print('Running save')
save(dge_results, file = 'dge_results.RData')

# F4
print('Running splot_diagnostic_plotsave')
futuriandgeDownstream::plot_diagnostic_plots(count_matrix_raw = counts_no_mirnas,
                      metadata_raw = sample_metadata)

# F5
print('futuriandgeDownstream::render_dge_html_report')
contrasts_for_report <- dplyr::mutate(contrasts, compared_gorups = paste0(
    studied_effect,
    '_vs_',
    baseline
  ))

for (groups in contrasts_for_report$compared_gorups) {

  sample_metadata$sample <- gsub("_T1",
                                 "",
                                 sample_metadata$sample_id)

  futuriandgeDownstream::render_dge_html_report(
    comparison_name = groups,
    dge_results_in = data.frame(dge_results$dge_results[[groups]]),
    metadata = sample_metadata,
    norm_counts = dge_results$normalised_counts,
    raw_counts = dge_results$raw_counts_no_gtf_miRNA,
    organism = organism,
    output_name = paste0(groups, ".html")
  )
}

