# Generates a small SYNTHETIC example dataset with the same sample design
# as the real study (6 ROIs, 72 samples, 2 control + 2 ethanol animals), so
# the pipeline can be run and explored before the real data is public.
# Counts are simulated, not real measurements, and gene IDs are placeholders
# (Gene0001, Gene0002, ...).
#
# Output is already committed under data/example/, so you don't need to run
# this to use the pipeline -- only to regenerate it.

suppressPackageStartupMessages(library(writexl))
set.seed(42)

# Sample design: mirrors the real study's ROI sizes and group splits.
roi_design <- list(
  list(roi = "plaqueCortex",    tissue = "brain", region = "cortex",      pathology = "plaque-bearing", n_ctrl = 6, n_etoh = 6),
  list(roi = "nonplaqueCortex", tissue = "brain", region = "cortex",      pathology = "plaque-free",    n_ctrl = 4, n_etoh = 4),
  list(roi = "plaqueHippo",     tissue = "brain", region = "hippocampus", pathology = "plaque-bearing", n_ctrl = 6, n_etoh = 6),
  list(roi = "nonplaqueHippo",  tissue = "brain", region = "hippocampus", pathology = "plaque-free",    n_ctrl = 4, n_etoh = 4),
  list(roi = "periportal",      tissue = "liver", region = "periportal",  pathology = NA,                n_ctrl = 8, n_etoh = 8),
  list(roi = "perivenous",      tissue = "liver", region = "perivenous",  pathology = NA,                n_ctrl = 8, n_etoh = 8)
)

# Roughly mirrors the manuscript's DEG pattern: liver >> hippocampus > cortex.
n_deg <- c(plaqueCortex = 0, nonplaqueCortex = 3, plaqueHippo = 40,
           nonplaqueHippo = 0, periportal = 80, perivenous = 300)

# 2 control + 2 ethanol "animals", each contributing to every ROI, same as
# the real study -- so example PCA plots can show batch (animal) alongside
# Control vs Alcohol-Treated. Animal is nested within treatment (each
# animal is entirely one group), matching the real design.
control_animals <- c("9001", "9002")
ethanol_animals <- c("9003", "9004")

samples <- do.call(rbind, lapply(roi_design, function(r) {
  treatment <- c(rep("control", r$n_ctrl), rep("ethanol", r$n_etoh))
  animal <- ifelse(treatment == "control",
                    rep_len(control_animals, r$n_ctrl)[cumsum(treatment == "control")],
                    rep_len(ethanol_animals, r$n_etoh)[cumsum(treatment == "ethanol")])
  data.frame(
    sample_id = sprintf("EX-%s-%s-%s-%02d", r$roi, treatment, animal, seq_along(treatment)),
    roi = r$roi, tissue = r$tissue, region = r$region, pathology = r$pathology,
    treatment = treatment, animal = animal, genotype = "example-genotype", sex = "male",
    stringsAsFactors = FALSE
  )
}))

# Simulate counts: negative-binomial background noise for every gene
# (size = 8, realistic overdispersion), plus a fold-change bump in each
# ROI's designated "DEG" genes with tighter dispersion (size = 40) so the
# injected signal is actually visible in a heatmap, not buried in noise.
n_genes <- 4000
gene_ids <- sprintf("Gene%04d", seq_len(n_genes))
gene_baseline <- rlnorm(n_genes, meanlog = 4, sdlog = 1.2)
bg_size <- 8
deg_size <- 40

counts <- matrix(NA_integer_, nrow = n_genes, ncol = nrow(samples),
                  dimnames = list(gene_ids, samples$sample_id))

for (r in names(n_deg)) {
  deg_genes <- sample(gene_ids, n_deg[[r]])
  deg_idx <- match(deg_genes, gene_ids)
  fold_change <- sample(c(-1, 1), n_deg[[r]], replace = TRUE) * runif(n_deg[[r]], 1.5, 3)

  roi_rows <- which(samples$roi == r)
  for (i in roi_rows) {
    mu <- gene_baseline
    gene_size <- rep(bg_size, n_genes)
    gene_size[deg_idx] <- deg_size
    if (samples$treatment[i] == "ethanol") {
      mu[deg_idx] <- mu[deg_idx] * 2^fold_change
    }
    counts[, samples$sample_id[i]] <- rnbinom(n_genes, mu = mu, size = gene_size)
  }
}

counts <- rbind(counts, "NegProbe-WTX" = rpois(ncol(counts), lambda = 2))
write.csv(counts, "data/example/RawCounts_example.csv")

# Metadata workbook, same shape 01_prepare_metadata.R already expects to
# parse: a few boilerplate rows, then a "*library name" header row, then
# one row per sample. Animal batch isn't a separate column here either --
# same as the real metadata, it's recovered from the sample_id in
# 01_prepare_metadata.R.
header <- c("*library name", "**tissue", "genotype", "treatment", "region", "pathology", "sex")
rows <- c(
  list("# Example dataset for pipeline demonstration"),
  list("# Synthetic data -- not the real study's measurements"),
  list(""),
  list(header),
  lapply(seq_len(nrow(samples)), function(i) {
    with(samples[i, ], c(sample_id, tissue, genotype, treatment, region, pathology, sex))
  })
)

width <- max(lengths(rows))
rows <- lapply(rows, function(r) { length(r) <- width; r })
meta_sheet <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)

write_xlsx(list(Metadata = meta_sheet), "data/example/MetaData_example.xlsx", col_names = FALSE)

cat(sprintf("Wrote data/example/RawCounts_example.csv (%d genes x %d samples)\n", nrow(counts), ncol(counts)))
cat("Wrote data/example/MetaData_example.xlsx\n")
