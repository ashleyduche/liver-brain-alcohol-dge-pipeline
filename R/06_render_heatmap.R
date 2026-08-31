# One heatmap panel per ROI (that ROI's own significant genes only, not a
# combined gene union across ROIs), arranged in a 2-column grid -- so 3
# included ROIs lay out as 2 panels on top, 1 below.

figures_dir <- file.path(cfg$paths$results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (length(included_rois) == 0) {
  cat("No eligible ROIs -- skipping heatmap.\n")
} else {
  # Collect each panel's data first, then bind once at the end -- binding
  # incrementally would silently drop the per-panel gene ordering below
  # (bind_rows can't preserve factor levels against an empty starting frame).
  roi_heatmap_pieces <- list()
  max_genes <- 0
  for (roi_name in included_rois) {
    sig <- roi_results[[roi_name]]$sig_results
    genes <- sig$gene[order(-sig$log2FoldChange)] # up-regulated genes together, then down -- stands in for clustering
    max_genes <- max(max_genes, length(genes))

    z <- t(scale(t(roi_results[[roi_name]]$norm_log2[genes, , drop = FALSE])))
    df <- as.data.frame(z)
    df$gene <- factor(rownames(z), levels = rev(genes))
    df <- pivot_longer(df, -gene, names_to = "sample_id", values_to = "z")
    df$ROI <- roi_results[[roi_name]]$roi_def$label

    roi_heatmap_pieces[[roi_name]] <- df
  }
  roi_heatmap_long <- bind_rows(roi_heatmap_pieces)

  p <- ggplot(roi_heatmap_long, aes(sample_id, gene, fill = z)) +
    geom_tile() +
    facet_wrap(~ ROI, scales = "free", ncol = 2) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "z-score") + # same red/blue as the volcano plot
    labs(title = "Significant DEGs by ROI (z-scored log2 normalized counts)", x = NULL, y = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.y = if (max_genes <= 60) element_text(size = 5) else element_blank())

  n_grid_rows <- ceiling(length(included_rois) / 2)
  ggsave(file.path(figures_dir, "significant_DEGs_heatmap_by_roi.pdf"), p, width = 9, height = 4.5 * n_grid_rows, limitsize = FALSE)
  ggsave(file.path(figures_dir, "significant_DEGs_heatmap_by_roi.png"), p, width = 9, height = 4.5 * n_grid_rows, dpi = 150, limitsize = FALSE)

  cat(sprintf("Wrote per-ROI heatmap (%d ROIs) to %s\n", length(included_rois), figures_dir))
}
