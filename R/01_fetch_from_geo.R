# Downloads the raw counts file from GEO (GSE324193) and replaces the local copy
# Only runs when data_source is geo in config.yml
#
# Metadata still comes from the local MetaData.xlsx rather than GEO
# GEO stores the sample details as free text linked to GSM accessions
# Those are not the *library name values used as column headers in the counts file

if (identical(cfg$data_source, "geo")) {
  suppressPackageStartupMessages(library(GEOquery))

  if (is.null(cfg$geo$accession) || is.na(cfg$geo$accession)) {
    stop("data_source is 'geo' but config.yml's geo.accession is not set. ",
         "Fill it in once GEO assigns an accession after publication.")
  }

  dir.create(dirname(cfg$paths$raw_counts_csv), recursive = TRUE, showWarnings = FALSE)

  # Download the supplementary files and find the counts file among them
  supp_files <- getGEOSuppFiles(cfg$geo$accession, baseDir = tempdir(), makeDirectory = TRUE)
  supp_path <- rownames(supp_files)[grepl(cfg$geo$supplementary_file, rownames(supp_files), fixed = TRUE)]
  stopifnot("Expected supplementary counts file not found in the GEO series" = length(supp_path) == 1)

  file.copy(supp_path, cfg$paths$raw_counts_csv, overwrite = TRUE)
  cat(sprintf("Fetched %s from GEO accession %s -> %s\n",
              cfg$geo$supplementary_file, cfg$geo$accession, cfg$paths$raw_counts_csv))
}
