# PCA per ROI on normalized counts, run before DESeq2 model fitting -- a QC
# check on the actual data going into the model, not a byproduct of it.
# Colored by Group, shaped by batch (animal). Animal is nested within Group
# here (each animal is entirely one treatment), so this checks for
# animal-driven separation, not a clean batch effect test.

figures_dir <- file.path(cfg$paths$results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

pca_plots <- list()

for (roi_name in names(roi_dds)) {
  norm_log2 <- roi_dds[[roi_name]]$norm_log2

  pca <- prcomp(t(norm_log2), scale. = FALSE)
  var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

  pca_df <- as.data.frame(pca$x[, 1:2])
  pca_df$sample_id <- rownames(pca_df)
  pca_df <- left_join(pca_df, meta[, c("sample_id", "Group", "batch")], by = "sample_id")

  p <- ggplot(pca_df, aes(PC1, PC2, color = Group, shape = batch)) +
    geom_point(size = 3) +
    labs(title = roi_dds[[roi_name]]$roi_def$label,
         x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         shape = "Batch (animal)") +
    theme_bw()

  pca_plots[[roi_name]] <- p
  ggsave(file.path(figures_dir, paste0(roi_name, "_PCA.pdf")), p, width = 6, height = 5)
  ggsave(file.path(figures_dir, paste0(roi_name, "_PCA.png")), p, width = 6, height = 5, dpi = 150)
}

cat(sprintf("Wrote %d PCA plots to %s\n", length(pca_plots), figures_dir))
