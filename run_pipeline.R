# Top-level driver for the liver & brain alcohol-treatment DGE pipeline.
# Run from the project root: Rscript run_pipeline.R

source("R/00_setup.R")
source("R/00_fetch_from_geo.R")             # no-op unless data_source: geo in config.yml
source("R/01_prepare_metadata.R")
source("R/02_run_dge_all_rois.R")           # loop-based; swap in 02_run_dge_all_rois_no_loop.R for the unlooped version
source("R/03_compare_deg_lists.R")
source("R/04_render_heatmap.R")
source("R/05_pca_plots.R")

cat("\nPipeline complete. See results/ for outputs:\n")
cat("  results/normalized_counts/  - per-ROI log2 DESeq2-normalized counts\n")
cat("  results/deg_tables/         - per-ROI full + significant DEG tables\n")
cat("  results/comparisons/        - shared/unique gene summaries across ROIs\n")
cat("  results/figures/            - PCA plots, significant-DEG heatmaps\n")
