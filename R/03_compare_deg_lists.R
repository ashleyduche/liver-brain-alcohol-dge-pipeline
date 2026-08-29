# Compares significant DEG lists across ROIs: shared (>=2 ROIs) vs unique
# (1 ROI). Only ROIs with >= min_sig_genes_downstream genes are included.

comparisons_dir <- file.path(cfg$paths$results_dir, "comparisons")
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

included_rois <- names(roi_results)[vapply(roi_results, `[[`, logical(1), "included")]
excluded_rois <- setdiff(names(roi_results), included_rois)

if (length(excluded_rois) > 0) {
  cat(sprintf("Excluding %s from cross-ROI comparison (< %d significant DEGs).\n",
              paste(excluded_rois, collapse = ", "), cfg$thresholds$min_sig_genes_downstream))
}

# One row per (gene, ROI) significant in that ROI, with tissue attached.
roi_tissue <- tibble(roi = names(cfg$rois),
                      tissue = vapply(cfg$rois, `[[`, character(1), "tissue"))

sig_long <- roi_results[included_rois] %>%
  lapply(function(r) r$sig_results["gene"]) %>%
  bind_rows(.id = "roi") %>%
  left_join(roi_tissue, by = "roi")

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
n_brain_rois <- length(intersect(included_rois, roi_tissue$roi[roi_tissue$tissue == "brain"]))
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
n_liver_rois <- length(intersect(included_rois, roi_tissue$roi[roi_tissue$tissue == "liver"]))
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
