# Builds each ROI's DESeqDataSet and estimates size factors -- normalized
# counts only, no differential testing yet. Splitting this out lets PCA/QC
# (04_pca_plots.R) run on normalized data before DESeq2 model fitting
# (05_run_dge_all_rois.R), which reuses the dds objects built here.

dir.create(cfg$paths$results_dir, recursive = TRUE, showWarnings = FALSE)
norm_dir <- file.path(cfg$paths$results_dir, "normalized_counts")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)

roi_dds <- list()

for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]
  roi_data <- subset_counts_by_tissue(meta, raw_counts, roi_def)
  coldata  <- roi_data$coldata
  counts   <- roi_data$counts

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts),
    colData   = coldata,
    design    = ~ Group
  )
  dds <- estimateSizeFactors(dds)
  norm_log2 <- log2(counts(dds, normalized = TRUE, replaced = FALSE) + 1)
  write.csv(norm_log2, file.path(norm_dir, paste0(roi_name, "_log2NormalizedCounts.csv")))

  roi_dds[[roi_name]] <- list(dds = dds, roi_def = roi_def, norm_log2 = norm_log2)
}

cat(sprintf("Normalized %d ROIs -> %s\n", length(roi_dds), norm_dir))
