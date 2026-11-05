# Parse the GEO metadata workbook into a clean table with one row per sample
# Adds the Control and Alcohol-Treated Group factor for DESeq2
# Also adds a batch column (animal ID) for PCA and QC
# Animal is not its own field in the metadata, but it is part of every sample_id (like Tg-Cont-6876-PlaqueHIP-01)
# so it is pulled out with a regex

raw_meta <- read_excel(cfg$paths$metadata_xlsx, sheet = cfg$paths$metadata_sheet, col_names = FALSE)

# Find the SAMPLES section header by name
header_row <- which(raw_meta[[1]] == "*library name")
stopifnot("Could not uniquely locate the '*library name' header row" = length(header_row) == 1)
colnames(raw_meta) <- as.character(raw_meta[header_row, ])

meta <- raw_meta %>%
  select(`*library name`, `**tissue`, genotype, treatment, region, pathology, sex) %>%
  dplyr::slice((header_row + 1):n()) %>%
  filter(!is.na(`**tissue`)) %>%
  transmute(
    sample_id = `*library name`,
    tissue    = `**tissue`,
    genotype,
    treatment,
    region,
    pathology,
    sex,
    batch = str_extract(sample_id, "[0-9]{4}"),
    Group = if_else(tolower(treatment) == tolower(cfg$group$treated_label),
                     cfg$group$treated_group, cfg$group$control_group),
    Group = factor(Group, levels = c(cfg$group$control_group, cfg$group$treated_group))
  )

raw_counts <- load_raw_counts(cfg$paths$raw_counts_csv, cfg$paths$negative_probe)

# Check that the metadata and counts contain the same samples
missing_in_counts <- setdiff(meta$sample_id, colnames(raw_counts))
missing_in_meta   <- setdiff(colnames(raw_counts), meta$sample_id)

if (length(missing_in_counts) > 0) {
  stop("Samples in metadata but missing from raw counts: ",
       paste(missing_in_counts, collapse = ", "))
}

if (length(missing_in_meta) > 0) {
  stop("Samples in raw counts but missing from metadata: ",
       paste(missing_in_meta, collapse = ", "))
}

# Save the parsed metadata
dir.create(dirname(cfg$paths$processed_metadata_csv), recursive = TRUE, showWarnings = FALSE)
write.csv(meta, cfg$paths$processed_metadata_csv, row.names = FALSE)

# capture.output turns the table into text so cat prints the summary line and table together
cat(sprintf("Parsed metadata for %d samples across %d ROIs.\n", nrow(meta), length(cfg$rois)),
    paste(capture.output(table(meta$region, meta$pathology, meta$Group, useNA = "ifany")), collapse = "\n"), "\n",
    sep = "")
