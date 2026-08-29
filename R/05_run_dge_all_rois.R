# Runs DESeq2 (Control vs Alcohol-Treated) for each ROI, reusing the
# already-normalized dds objects from 03_normalize_counts.R. Same steps as
# the original per-ROI scripts, looped instead of copy-pasted 6 times.
# Prefer no loop? See 05_run_dge_all_rois_no_loop.R.

deg_dir <- file.path(cfg$paths$results_dir, "deg_tables")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

roi_results <- list()

for (roi_name in names(roi_dds)) {
  roi_def <- roi_dds[[roi_name]]$roi_def
  dds <- roi_dds[[roi_name]]$dds
  cat(sprintf("\n=== %s (%s) ===\n", roi_name, roi_def$label))
  cat(sprintf("  %d samples (%d Control / %d Alcohol-Treated)\n",
              ncol(dds),
              sum(dds$Group == cfg$group$control_group),
              sum(dds$Group == cfg$group$treated_group)))

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
    norm_log2   = roi_dds[[roi_name]]$norm_log2,
    n_sig       = n_sig,
    included    = included
  )
}

n_excluded <- sum(!vapply(roi_results, `[[`, logical(1), "included"))
cat(sprintf("\n%d/%d ROIs pass the >= %d significant-gene threshold for downstream comparison.\n",
            length(roi_results) - n_excluded, length(roi_results),
            cfg$thresholds$min_sig_genes_downstream))
