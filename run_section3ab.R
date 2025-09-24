# Safe runner for sections 3a-b: read data, run DESeq2, save results and volcano
options(stringsAsFactors = FALSE)

pkg_check <- function(pkgs) {
  to_install <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(to_install)) {
    message("Installing missing packages: ", paste(to_install, collapse=", "))
    install.packages(to_install, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}

bioc_check <- function(pkgs) {
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) {
    if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
    BiocManager::install(need, update = FALSE)
  }
}

pkg_check(c("dplyr","readr","stringr","magrittr","tibble","ggplot2","ggpubr","jsonlite"))
bioc_check(c("DESeq2","apeglm","EnhancedVolcano"))

# load packages
library(readr); library(dplyr); library(stringr); library(magrittr); library(tibble)
library(DESeq2); library(apeglm); library(EnhancedVolcano); library(ggplot2)

# paths
metadata_path <- file.path("data","SRP192714","metadata_SRP192714.tsv")
expr_path     <- file.path("data","SRP192714","SRP192714.tsv")
results_dir <- "results"
plots_dir <- "plots"
if (!dir.exists(results_dir)) dir.create(results_dir)
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# helper for clean stop
err <- function(msg) { stop(msg, call. = FALSE) }

# Read files
if (!file.exists(metadata_path)) err(paste("Missing metadata file:", metadata_path))
if (!file.exists(expr_path)) err(paste("Missing expression file:", expr_path))

metadata <- readr::read_tsv(metadata_path, show_col_types = FALSE)
expression_df <- readr::read_tsv(expr_path, show_col_types = FALSE)

# Parse command line args: support --ignore-mutation to skip using mutation_status in design
args <- commandArgs(trailingOnly = TRUE)
ignore_mutation <- any(grepl("--ignore-mutation", args))
if (ignore_mutation) message('Running with --ignore-mutation: design will not use mutation_status')

# Ensure Gene column exists
if (!"Gene" %in% colnames(expression_df)) err("Expression file must have a column named 'Gene'")

# Make rownames
expression_df <- expression_df %>% tibble::column_to_rownames(var = "Gene")

# Check sample names in metadata
if (!"refinebio_accession_code" %in% colnames(metadata)) {
  err("metadata is missing column 'refinebio_accession_code'")
}

# Check for exact matches
missing_cols <- setdiff(metadata$refinebio_accession_code, colnames(expression_df))
extra_cols   <- setdiff(colnames(expression_df), metadata$refinebio_accession_code)
if (length(missing_cols) > 0) {
  message('Samples in metadata not found in expression columns:')
  print(missing_cols)
  err('Column name mismatch: fix sample names in metadata or expression file')
}
# Reorder expression columns to metadata order
expression_df <- expression_df %>% dplyr::select(all_of(metadata$refinebio_accession_code))

# Prepare mutation_status unless user requested to ignore it
labels_path <- file.path("sample_labels.csv")
if (!ignore_mutation) {
  if (!"refinebio_title" %in% colnames(metadata)) {
    message("metadata missing 'refinebio_title' column; mutation_status will be NA")
    metadata$refinebio_title <- NA_character_
  }
  # automatic detection from titles
  metadata <- metadata %>%
    dplyr::mutate(mutation_status = dplyr::case_when(
      stringr::str_detect(refinebio_title, regex("R98S", ignore_case = TRUE)) ~ "R98S",
      stringr::str_detect(refinebio_title, regex("\\b(WT|reference)\\b", ignore_case = TRUE)) ~ "reference",
      TRUE ~ NA_character_
    ))

  # Allow explicit mapping file to override detection
  if (file.exists(labels_path)) {
    message('Reading sample labels from ', labels_path)
    labels <- readr::read_csv(labels_path, show_col_types = FALSE)
    if (!all(c('refinebio_accession_code','mutation_status') %in% names(labels))) {
      stop('sample_labels.csv must contain columns: refinebio_accession_code, mutation_status')
    }
    metadata <- metadata %>% dplyr::left_join(labels, by = 'refinebio_accession_code') %>%
      dplyr::mutate(mutation_status = dplyr::coalesce(mutation_status.y, mutation_status.x)) %>%
      dplyr::select(-mutation_status.x, -mutation_status.y)
  } else {
    # create a template for the user to fill if any mutation_status are NA after detection
    if (any(is.na(metadata$mutation_status))) {
      template <- metadata %>% dplyr::select(refinebio_accession_code) %>%
        dplyr::mutate(mutation_status = NA_character_)
      readr::write_csv(template, labels_path)
      stop(paste0('No sample_labels.csv found. A template has been written to ', labels_path,
                  '. Please fill the mutation_status values (e.g. reference or R98S) and re-run the script.'))
    }
  }

  # ensure factor levels
  metadata$mutation_status <- factor(metadata$mutation_status, levels = c("reference","R98S"))
} else {
  # When ignoring mutation_status, ensure the column exists (can be NA) but will not be used in design
  if (!"mutation_status" %in% colnames(metadata)) metadata$mutation_status <- NA_character_
}

# Filter low counts
filtered_expression_df <- expression_df[rowSums(as.matrix(expression_df), na.rm = TRUE) >= 10, , drop = FALSE]
if (nrow(filtered_expression_df) == 0) err('No genes passed the count filter (rowSums >=10)')

# Convert to integer matrix
gene_matrix <- round(as.matrix(filtered_expression_df))

# Reorder metadata rows to match columns
metadata <- metadata[match(colnames(gene_matrix), metadata$refinebio_accession_code), , drop = FALSE]
if (!all.equal(as.character(metadata$refinebio_accession_code), colnames(gene_matrix))) err('After reordering, sample names still do not match')

# Build DESeq2 object
design_formula <- if (ignore_mutation) as.formula('~ 1') else as.formula('~ mutation_status')
dds <- DESeqDataSetFromMatrix(countData = gene_matrix, colData = metadata, design = design_formula)

# make warnings print immediately
old_warn <- options("warn")
options(warn = 1)
if (requireNamespace("BiocParallel", quietly = TRUE)) {
  BiocParallel::register(BiocParallel::SerialParam())
}

# Run DESeq with warning capture
dds <- tryCatch({
  withCallingHandlers({
    DESeq(dds)
  }, warning = function(w) {
    message('DESeq warning: ', conditionMessage(w))
    invokeRestart('muffleWarning')
  })
}, error = function(e) { options(warn = old_warn[[1]]); err(paste('DESeq failed:', e$message)) })
message('DESeq finished — continuing script')
options(warn = old_warn[[1]])

if (!ignore_mutation) {
  res <- results(dds)
  # Use apeglm shrink
  res_shrunk <- tryCatch({ lfcShrink(dds, coef = 2, type = "apeglm", res = res) }, error = function(e) {
    message('apeglm lfcShrink failed, falling back to results without shrink: ', e$message)
    res
  })

  # Make deseq_df
  deseq_df <- as.data.frame(res_shrunk) %>% tibble::rownames_to_column('Gene') %>%
    dplyr::mutate(threshold = ifelse(is.na(padj), FALSE, padj < 0.05)) %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))

  # Save TSV
  out_tsv <- file.path(results_dir, "SRP192714_diff_expr_results.tsv")
  readr::write_tsv(deseq_df, out_tsv)
  message('Wrote ', out_tsv)

  # Create volcano plot
  volcano_plot <- EnhancedVolcano::EnhancedVolcano(
    deseq_df,
    lab = deseq_df$Gene,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.01
  )

  # Save volcano
  out_png <- file.path(plots_dir, 'SRP192714_volcano_plot.png')
  tryCatch({ ggplot2::ggsave(plot = volcano_plot, filename = out_png, width = 8, height = 6) ; message('Wrote ', out_png) }, error = function(e) { message('Failed to save volcano PNG: ', e$message) })

  # Save objects to RDS for inspection
  saveRDS(deseq_df, file = file.path(results_dir, 'deseq_df.rds'))
  saveRDS(volcano_plot, file = file.path(results_dir, 'volcano_plot.rds'))

} else {
  # When ignoring mutation_status, write normalized counts and save dds object for downstream analysis
  vsd <- tryCatch({ vst(dds) }, error = function(e) { message('vst failed, trying rlog: ', e$message); tryCatch({ rlog(dds) }, error = function(e2) { err(paste('vst/rlog both failed:', e2$message)) }) })
  norm_mat <- assay(vsd)
  out_norm <- file.path(results_dir, 'SRP192714_normalized_counts.tsv')
  readr::write_tsv(as.data.frame(norm_mat) %>% tibble::rownames_to_column('Gene'), out_norm)
  message('Wrote normalized counts to ', out_norm)
  # save the DESeqDataSet
  saveRDS(dds, file = file.path(results_dir, 'SRP192714_dds_ignore_mutation.rds'))
  message('Saved DESeqDataSet to results/SRP192714_dds_ignore_mutation.rds')
}

message('Script completed successfully')
