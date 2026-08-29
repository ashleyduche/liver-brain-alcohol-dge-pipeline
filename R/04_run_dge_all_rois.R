# Runs DESeq2 (Control vs Alcohol-Treated) for each ROI in config.yml.
# Same steps as the original per-ROI scripts, looped instead of
# copy-pasted 6 times. Prefer no loop? See 04_run_dge_all_rois_no_loop.R.

dir.create(cfg$paths$results_dir, recursive = TRUE, showWarnings = FALSE)
norm_dir <- file.path(cfg$paths$results_dir, "normalized_counts")
deg_dir  <- file.path(cfg$paths$results_dir, "deg_tables")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

roi_results <- list()

for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]
  cat(sprintf("\n=== %s (%s) ===\n", roi_name, roi_def$label))

  roi_data <- subset_counts_by_tissue(meta, raw_counts, roi_def)
  coldata  <- roi_data$coldata
  counts   <- roi_data$counts

  cat(sprintf("  %d samples (%d Control / %d Alcohol-Treated)\n",
              nrow(coldata),
              sum(coldata$Group == cfg$group$control_group),
              sum(coldata$Group == cfg$group$treated_group)))

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(counts),
    colData   = coldata,
    design    = ~ Group
  )

  dds <- estimateSizeFactors(dds)
  norm_log2 <- log2(counts(dds, normalized = TRUE, replaced = FALSE) + 1)
  write.csv(norm_log2, file.path(norm_dir, paste0(roi_name, "_log2NormalizedCounts.csv")))

  keep_genes <- rowSums(counts(dds)) >= cfg$thresholds$min_rowsum_prefilter
  dds <- dds[keep_genes, ]

  dds <- DESeq(dds)
  res <- results(dds)

  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  write.csv(res_df, file.path(deg_dir, paste0(roi_name, "_AllDEGs.csv")), row.names = FALSE)

  res_df_no_na <- res_df[!is.na(res_df$padj), ] # DESeq2's own independent filtering already applied

  sig_df <- res_df_no_na %>%
    filter(log2FoldChange > cfg$thresholds$abs_log2fc | log2FoldChange < -cfg$thresholds$abs_log2fc,
           padj <= cfg$thresholds$padj)
  write.csv(sig_df, file.path(deg_dir, paste0(roi_name, "_SigDEGs_padj", cfg$thresholds$padj,
                                                "_absLFC", cfg$thresholds$abs_log2fc, ".csv")),
            row.names = FALSE)

  n_sig <- nrow(sig_df)
  included <- n_sig >= cfg$thresholds$min_sig_genes_downstream
  cat(sprintf("  %d significant DEGs (padj <= %s, |log2FC| > %s) -> %s\n",
              n_sig, cfg$thresholds$padj, cfg$thresholds$abs_log2fc,
              if (included) "included in cross-ROI comparison" else "EXCLUDED (< min gene threshold)"))

  roi_results[[roi_name]] <- list(
    roi_def     = roi_def,
    sig_results = sig_df,
    norm_log2   = norm_log2,
    n_sig       = n_sig,
    included    = included
  )
}

n_excluded <- sum(!vapply(roi_results, `[[`, logical(1), "included"))
cat(sprintf("\n%d/%d ROIs pass the >= %d significant-gene threshold for downstream comparison.\n",
            length(roi_results) - n_excluded, length(roi_results),
            cfg$thresholds$min_sig_genes_downstream))
