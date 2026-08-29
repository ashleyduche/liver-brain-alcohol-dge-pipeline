# No-loop alternative to R/02_run_dge_all_rois.R -- same DESeq2 steps for
# ONE ROI. The steps don't change by tissue, so just interchange roi_name
# below and re-run to analyze a different region.

roi_name <- "plaqueCortex" # <- change to: nonplaqueCortex, plaqueHippo, nonplaqueHippo, periportal, perivenous

if (!exists("cfg"))  source("R/00_setup.R")
if (!exists("meta")) source("R/01_prepare_metadata.R")

dir.create(cfg$paths$results_dir, recursive = TRUE, showWarnings = FALSE)
norm_dir <- file.path(cfg$paths$results_dir, "normalized_counts")
deg_dir  <- file.path(cfg$paths$results_dir, "deg_tables")
dir.create(norm_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
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

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts), colData = coldata, design = ~ Group)

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
