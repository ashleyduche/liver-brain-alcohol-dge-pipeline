# Top level driver for the liver and brain alcohol treatment DGE pipeline
# Run from the project root with: Rscript run_pipeline.R

source("R/00_setup.R")
source("R/01_fetch_from_geo.R")             # only runs when data_source is geo in config.yml
source("R/02_prepare_metadata.R")
source("R/03_pca_plots.R")                  # raw count PCA and QC, before any DESeq2 object exists
source("R/04_run_dge_all_rois.R")           # loop version, or swap in 04_run_dge_all_rois_no_loop.R for one ROI at a time
source("R/05_compare_deg_lists.R")
source("R/06_render_heatmap.R")

cat("\nPipeline complete. See results/ for outputs:\n")
cat("  results/normalized_counts/  - per-ROI log2 DESeq2-normalized counts\n")
cat("  results/deg_tables/         - per-ROI full + significant DEG tables\n")
cat("  results/comparisons/        - shared/unique gene summaries across ROIs\n")
cat("  results/figures/            - PCA plots, significant-DEG heatmaps\n")
