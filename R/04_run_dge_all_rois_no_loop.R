# No-loop alternative to R/04_run_dge_all_rois.R -- same steps for ONE
# ROI, in the same order as the original per-ROI scripts: raw-count PCA,
# then build dds -> normalize -> filter -> DESeq2. The steps don't change
# by tissue, so just interchange roi_name below and re-run.

roi_name <- "plaqueCortex" # <- change to: nonplaqueCortex, plaqueHippo, nonplaqueHippo, periportal, perivenous

if (!exists("cfg"))  source("R/00_setup.R")
if (!exists("meta")) source("R/02_prepare_metadata.R")

dir.create(cfg$paths$results_dir, recursive = TRUE, showWarnings = FALSE)
norm_dir    <- file.path(cfg$paths$results_dir, "normalized_counts")
deg_dir     <- file.path(cfg$paths$results_dir, "deg_tables")
figures_dir <- file.path(cfg$paths$results_dir, "figures")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
if (!exists("roi_results")) roi_results <- list()

roi_def <- cfg$rois[[roi_name]]
if (is.null(roi_def)) {
  stop("Unknown roi_name '", roi_name, "'. Must be one of: ", paste(names(cfg$rois), collapse = ", "))
}

cat(sprintf("\n=== %s (%s) ===\n", roi_name, roi_def$label))

roi_data <- subset_counts_by_tissue(meta, raw_counts, roi_def)
coldata  <- roi_data$coldata
counts   <- roi_data$counts
cat(sprintf("  %d samples (%d Control / %d Alcohol-Treated)\n",
            nrow(coldata),
            sum(coldata$Group == cfg$group$control_group),
            sum(coldata$Group == cfg$group$treated_group)))

# PCA on raw counts, before DESeq2 -- same as the original per-ROI scripts.
# coldata is already in the same row order as counts, so Group/batch are
# just pulled straight off it -- no join needed.
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

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = ~ Group)

dds <- estimateSizeFactors(dds)
norm_log2 <- log2(counts(dds, normalized = TRUE, replaced = FALSE) + 1)
write.csv(norm_log2, file.path(norm_dir, paste0(roi_name, "_log2NormalizedCounts.csv")))

keep_genes <- rowSums(counts(dds)) >= cfg$thresholds$min_rowsum_prefilter
dds <- dds[keep_genes, ]

dds <- DESeq(dds)
res <- results(dds)

# MA plot -- same as the original per-ROI scripts (plotMA(res)).
pdf(file.path(figures_dir, paste0(roi_name, "_MA.pdf")))
plotMA(res, main = roi_def$label)
dev.off()
png(file.path(figures_dir, paste0(roi_name, "_MA.png")), width = 800, height = 700, res = 120)
plotMA(res, main = roi_def$label)
dev.off()

# Volcano plot: log2FoldChange vs -log10(padj), significant genes in red.
is_sig <- !is.na(res$padj) & res$padj <= cfg$thresholds$padj &
  (res$log2FoldChange > cfg$thresholds$abs_log2fc | res$log2FoldChange < -cfg$thresholds$abs_log2fc)
volcano_col <- ifelse(is_sig, "red", "grey60")

pdf(file.path(figures_dir, paste0(roi_name, "_Volcano.pdf")))
plot(res$log2FoldChange, -log10(res$padj), col = volcano_col, pch = 20,
     xlab = "log2 Fold Change", ylab = "-log10(padj)", main = paste0(roi_def$label, ": Volcano Plot"))
abline(v = c(-cfg$thresholds$abs_log2fc, cfg$thresholds$abs_log2fc), lty = 2, col = "grey40")
abline(h = -log10(cfg$thresholds$padj), lty = 2, col = "grey40")
dev.off()
png(file.path(figures_dir, paste0(roi_name, "_Volcano.png")), width = 800, height = 700, res = 120)
plot(res$log2FoldChange, -log10(res$padj), col = volcano_col, pch = 20,
     xlab = "log2 Fold Change", ylab = "-log10(padj)", main = paste0(roi_def$label, ": Volcano Plot"))
abline(v = c(-cfg$thresholds$abs_log2fc, cfg$thresholds$abs_log2fc), lty = 2, col = "grey40")
abline(h = -log10(cfg$thresholds$padj), lty = 2, col = "grey40")
dev.off()

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
