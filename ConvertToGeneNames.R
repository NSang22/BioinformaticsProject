# This script is not used any more

print("Running ConvertToGeneNames.R")

# uncomment the below to download package dependencies
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager", repos = "https://cloud.r-project.org")
# BiocManager::install(version = "3.21")
# install.packages("tidyverse", repos = "https://cloud.r-project.org")
# install.packages("devtools", repos = "https://cloud.r-project.org")
# install.packages("GenomeInfoDb", repos = "https://cloud.r-project.org")
# Install the Homo Sapiens package
# if (!("org.Hs.eg.db" %in% installed.packages())) {
#   # Install this package if it isn't installed yet
#   BiocManager::install("org.Hs.eg.db", update = FALSE)
# }
# # Attach the library
# library(org.Hs.eg.db)
# # We will need this so we can use the pipe: %>%
# library(magrittr)

# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

# Define the file path to the data directory
data_dir <- file.path("data", "SRP192714")

# Declare the file path to the gene expression matrix file
# inside directory saved as `data_dir`
data_file <- file.path(data_dir, "SRP192714.tsv")

# Declare the file path to the metadata file
# inside the directory saved as `data_dir`
metadata_file <- file.path(data_dir, "metadata_SRP192714.tsv")



# Read in metadata TSV file
metadata <- readr::read_tsv(metadata_file)

# Read in data TSV file
expression_df <- readr::read_tsv(data_file) %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated Symbols
mapped_list <- mapIds(
  org.Hs.eg.db,
  keys = expression_df$Gene,
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "list"
)

# Let's use the `head()` function for a preview of our mapped list
head(mapped_list)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)

head(mapped_df)

# Use the `summary()` function to show the distribution of Symbol values
# We need to use `as.factor()` here to get the count of unique values
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Symbol), maxsum = 10)

multi_mapped <- mapped_df %>%
  # Let's count the number of times each Ensembl ID appears in `Ensembl` column
  dplyr::count(Ensembl, name = "symbol_count") %>%
  # Arrange by the genes with the highest number of Symbols mapped
  dplyr::arrange(desc(symbol_count))

# Let's look at the first 6 rows of our `multi_mapped` object
head(multi_mapped)

collapsed_mapped_df <- mapped_df %>%
  # Group by Ensembl IDs
  dplyr::group_by(Ensembl) %>%
  # Collapse the Symbol IDs `mapped_df` into one column named `all_symbol_ids`
  dplyr::summarize(all_symbols = paste(Symbol, collapse = ";"))

collapsed_mapped_df %>%
  # Filter `collapsed_mapped_df` to include only the rows where
  # `all_symbols` values include the ";" character --
  # these are the rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbols, ";")) %>%
  # We only need a preview here
  head()

final_mapped_df <- data.frame(
  "first_mapped_symbol" = mapIds(
    org.Hs.eg.db,
    keys = expression_df$Gene,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first" # Keep only the first mapped value for each Ensembl ID
  )
) %>%
  # Make an `Ensembl` column to store the rownames
  tibble::rownames_to_column("Ensembl") %>%
  # Add the multiple mappings data from `collapsed_mapped_df` using Ensembl IDs
  dplyr::inner_join(collapsed_mapped_df, by = "Ensembl") %>%
  # Now let's add on the rest of the expression data
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene"))

final_mapped_df %>%
  # Filter `final_mapped_df` to rows with multiple mapped values
  dplyr::filter(stringr::str_detect(all_symbols, ";")) %>%
  head()

# Write mapped and annotated data frame to output file
readr::write_tsv(final_mapped_df, file.path(
  results_dir,
  "SRP192714_Symbols.tsv"
))
