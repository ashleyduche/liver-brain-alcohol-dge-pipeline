# PCA per ROI, before any DESeq2 object is built -- same base R approach as
# the original per-ROI scripts: pca <- prcomp(t(log2(counts+1))), then
# plot(pca) for a scree plot and plot(pca$x[,1], pca$x[,2], col=...) colored
# Control (green) vs Alcohol-Treated (red). coldata is already in the same
# row order as counts (see subset_counts_by_tissue), so Group/batch are
# just pulled straight off coldata -- no join needed.

figures_dir <- file.path(cfg$paths$results_dir, "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

for (roi_name in names(cfg$rois)) {
  roi_def <- cfg$rois[[roi_name]]
  roi_data <- subset_counts_by_tissue(meta, raw_counts, roi_def)
  coldata  <- roi_data$coldata
  counts   <- roi_data$counts

  pca <- prcomp(t(log2(counts + 1)))

  point_col <- ifelse(coldata$Group == cfg$group$control_group, "green3", "red")
  batch_levels <- sort(unique(coldata$batch))
  point_pch <- match(coldata$batch, batch_levels)
  y_range <- range(pca$x[, 2])
  y_lim <- c(y_range[1], y_range[2] + diff(y_range) * 0.45) # headroom for the legend

  pdf(file.path(figures_dir, paste0(roi_name, "_PCA.pdf")))
  plot(pca)
  plot(pca$x[, 1], pca$x[, 2], ylim = y_lim,
       col = point_col, pch = point_pch, cex = 1.5,
       xlab = "PC1", ylab = "PC2",
       main = paste0(roi_def$label, ": Control (green) vs Alcohol-Treated (red)"))
  legend("topright", legend = batch_levels, pch = seq_along(batch_levels), title = "Batch (animal)", bty = "n")
  dev.off()

  png(file.path(figures_dir, paste0(roi_name, "_PCA.png")), width = 800, height = 700, res = 120)
  plot(pca$x[, 1], pca$x[, 2], ylim = y_lim,
       col = point_col, pch = point_pch, cex = 1.5,
       xlab = "PC1", ylab = "PC2",
       main = paste0(roi_def$label, ": Control (green) vs Alcohol-Treated (red)"))
  legend("topright", legend = batch_levels, pch = seq_along(batch_levels), title = "Batch (animal)", bty = "n")
  dev.off()
}

cat(sprintf("Wrote %d PCA plots to %s\n", length(cfg$rois), figures_dir))
