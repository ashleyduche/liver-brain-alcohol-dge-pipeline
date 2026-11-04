suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(DESeq2)
})

# Load the counts matrix as integers with gene symbols as rownames
# Kept as a function since DESeq2 needs a real matrix, not tidy data
load_raw_counts <- function(csv_path, negative_probe) {
  counts <- read.csv(csv_path, row.names = 1, check.names = FALSE)
  counts <- counts[rownames(counts) != negative_probe, , drop = FALSE]

  gene_names <- rownames(counts)
  counts <- apply(counts, 2, as.integer)
  rownames(counts) <- gene_names

  counts
}

# Subset metadata and counts to one ROI (tissue, region and pathology) with count columns in the same order
# Only data wrangling here, no DESeq2
subset_counts_by_tissue <- function(meta, counts, roi_def) {
  keep <- meta$tissue == roi_def$tissue & meta$region == roi_def$region
  if (!is.null(roi_def$pathology) && !is.na(roi_def$pathology)) {
    keep <- keep & meta$pathology == roi_def$pathology
  }

  coldata <- meta[keep, , drop = FALSE]
  roi_counts <- counts[, match(coldata$sample_id, colnames(counts)), drop = FALSE]
  stopifnot(identical(colnames(roi_counts), coldata$sample_id))

  list(counts = roi_counts, coldata = coldata)
}
