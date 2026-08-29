# Compares significant DEGs across the three ROIs that clear the >=10-gene
# threshold: plaqueHippo, periportal, perivenous (cortex ROIs never have
# enough DEGs to compare -- see the DESeq2 step above).

comparisons_dir <- file.path(cfg$paths$results_dir, "comparisons")
dir.create(comparisons_dir, recursive = TRUE, showWarnings = FALSE)

included_rois <- c("plaqueHippo", "periportal", "perivenous")

hip    <- roi_results$plaqueHippo$sig_results$gene
portal <- roi_results$periportal$sig_results$gene
venous <- roi_results$perivenous$sig_results$gene

cat("\n=== Cross-ROI DEG comparison ===\n")
cat(sprintf("hip: %d genes | portal: %d genes | venous: %d genes\n",
            length(hip), length(portal), length(venous)))

hip_portal    <- intersect(hip, portal)
hip_venous    <- intersect(hip, venous)
portal_venous <- intersect(portal, venous)
all_three     <- intersect(intersect(hip, portal), venous)

cat(sprintf("  hip & portal shared: %d\n", length(hip_portal)))
cat(sprintf("  hip & venous shared: %d\n", length(hip_venous)))
cat(sprintf("  portal & venous shared: %d\n", length(portal_venous)))
cat(sprintf("  all three shared: %d\n", length(all_three)))

# One table: every significant gene, and which ROI(s) it showed up in.
all_genes <- unique(c(hip, portal, venous))
gene_table <- data.frame(
  gene   = all_genes,
  hip    = all_genes %in% hip,
  portal = all_genes %in% portal,
  venous = all_genes %in% venous
)
gene_table$n_rois <- rowSums(gene_table[, c("hip", "portal", "venous")])
gene_table$status <- ifelse(gene_table$n_rois > 1, "shared", "unique")

write.csv(gene_table, file.path(comparisons_dir, "gene_summary.csv"), row.names = FALSE)
cat(sprintf("Total significant genes across all 3 ROIs: %d (%d shared, %d unique)\n",
            nrow(gene_table), sum(gene_table$status == "shared"), sum(gene_table$status == "unique")))
