# PCA per ROI on raw counts, before any DESeq2 object is built -- same
# approach as the original per-ROI scripts: pca <- prcomp(t(log2(counts+1))).
# Colored by Group (Control/Alcohol-Treated), shaped by batch (animal).

figures_dir <- file.path(cfg$paths$results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

pca_plots <- list()

for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]
  roi_data <- subset_counts_by_tissue(meta, raw_counts, roi_def)
  coldata  <- roi_data$coldata
  counts   <- roi_data$counts

  pca <- prcomp(t(log2(counts + 1)))
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, coldata[, c("sample_id", "Group", "batch")], by = "sample_id")

  p <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = batch)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(Control = "green3", "Alcohol-Treated" = "red")) +
    labs(title = roi_def$label,
         x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         shape = "Batch (animal)") +
    theme_bw()

  pca_plots[[roi_name]] <- p
  ggsave(file.path(figures_dir, paste0(roi_name, "_PCA.pdf")), p, width = 6, height = 5)
  ggsave(file.path(figures_dir, paste0(roi_name, "_PCA.png")), p, width = 6, height = 5, dpi = 150)
}

cat(sprintf("Wrote %d PCA plots to %s\n", length(pca_plots), figures_dir))
