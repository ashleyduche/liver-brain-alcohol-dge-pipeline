# Compares significant DEG lists across ROIs: shared (>=2 ROIs) vs unique
# (1 ROI). Only ROIs with >= min_sig_genes_downstream genes are included.

comparisons_dir <- file.path(cfg$paths$results_dir, "comparisons")
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

# Which ROIs have enough significant genes to compare across?
included_rois <- c()
excluded_rois <- c()
for (roi_name in names(roi_results)) {
  if (roi_results[[roi_name]]$included) {
    included_rois <- c(included_rois, roi_name)
  } else {
    excluded_rois <- c(excluded_rois, roi_name)
  }
}
if (length(excluded_rois) > 0) {
  cat(sprintf("Excluding %s from cross-ROI comparison (< %d significant DEGs).\n",
              paste(excluded_rois, collapse = ", "), cfg$thresholds$min_sig_genes_downstream))
}

# One row per (gene, ROI) that was significant in that ROI, with tissue attached.
sig_long <- tibble(gene = character(), roi = character(), tissue = character())
for (roi_name in included_rois) {
  genes  <- roi_results[[roi_name]]$sig_results$gene
  tissue <- cfg$rois[[roi_name]]$tissue
  sig_long <- bind_rows(sig_long, tibble(gene = genes, roi = roi_name, tissue = tissue))
}

cat("\n=== Cross-ROI DEG comparison ===\n")

# --- All eligible ROIs together (brain + liver) ---
gene_summary_all <- sig_long %>%
  distinct(gene, roi) %>%
  group_by(gene) %>%
  summarise(n_rois_significant = n(),
            rois = paste(sort(roi), collapse = ";"),
            status = if_else(n_rois_significant > 1, "shared", "unique"),
            .groups = "drop") %>%
  arrange(desc(n_rois_significant), gene)
write.csv(gene_summary_all, file.path(comparisons_dir, "all_rois_gene_summary.csv"), row.names = FALSE)
cat(sprintf("  [all ROIs] %d ROIs, %d total significant genes (%d shared, %d unique)\n",
            length(included_rois), nrow(gene_summary_all),
            sum(gene_summary_all$status == "shared"), sum(gene_summary_all$status == "unique")))

# --- Brain-only comparison ---
gene_summary_brain <- sig_long %>%
  filter(tissue == "brain") %>%
  distinct(gene, roi) %>%
  group_by(gene) %>%
  summarise(n_rois_significant = n(),
            rois = paste(sort(roi), collapse = ";"),
            status = if_else(n_rois_significant > 1, "shared", "unique"),
            .groups = "drop") %>%
  arrange(desc(n_rois_significant), gene)
write.csv(gene_summary_brain, file.path(comparisons_dir, "brain_rois_gene_summary.csv"), row.names = FALSE)

n_brain_rois <- 0
for (roi_name in included_rois) {
  if (cfg$rois[[roi_name]]$tissue == "brain") n_brain_rois <- n_brain_rois + 1
}
cat(sprintf("  [brain ROIs] %d ROIs, %d total significant genes (%d shared, %d unique)\n",
            n_brain_rois, nrow(gene_summary_brain),
            sum(gene_summary_brain$status == "shared"), sum(gene_summary_brain$status == "unique")))

# --- Liver-only comparison ---
gene_summary_liver <- sig_long %>%
  filter(tissue == "liver") %>%
  distinct(gene, roi) %>%
  group_by(gene) %>%
  summarise(n_rois_significant = n(),
            rois = paste(sort(roi), collapse = ";"),
            status = if_else(n_rois_significant > 1, "shared", "unique"),
            .groups = "drop") %>%
  arrange(desc(n_rois_significant), gene)
write.csv(gene_summary_liver, file.path(comparisons_dir, "liver_rois_gene_summary.csv"), row.names = FALSE)

n_liver_rois <- 0
for (roi_name in included_rois) {
  if (cfg$rois[[roi_name]]$tissue == "liver") n_liver_rois <- n_liver_rois + 1
}
cat(sprintf("  [liver ROIs] %d ROIs, %d total significant genes (%d shared, %d unique)\n",
            n_liver_rois, nrow(gene_summary_liver),
            sum(gene_summary_liver$status == "shared"), sum(gene_summary_liver$status == "unique")))

# --- Genes shared between brain and liver tissue-level DEG sets ---
if (n_brain_rois > 0 && n_liver_rois > 0) {
  tissue_summary <- sig_long %>%
    distinct(gene, tissue) %>%
    mutate(present = TRUE) %>%
    pivot_wider(names_from = tissue, values_from = present, values_fill = FALSE) %>%
    mutate(status = case_when(
      brain & liver ~ "shared_brain_liver",
      brain          ~ "brain_only",
      liver          ~ "liver_only"
    ))
  write.csv(tissue_summary, file.path(comparisons_dir, "brain_vs_liver_gene_summary.csv"), row.names = FALSE)
  cat(sprintf("  [brain vs liver] %d genes shared, %d brain-only, %d liver-only\n",
              sum(tissue_summary$status == "shared_brain_liver"),
              sum(tissue_summary$status == "brain_only"),
              sum(tissue_summary$status == "liver_only")))
}
