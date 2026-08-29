# Downloads the raw counts file from GEO once published, overwriting the
# local copy. No-op unless data_source: geo in config.yml.
#
# Metadata still comes from the local MetaData.xlsx, not GEO -- GEO only
# exposes it as free-text characteristics keyed to GSM accession/title, not
# to the *library name values the counts matrix uses as column headers, and
# that mapping isn't reliable to re-derive until a real accession exists.

if (identical(cfg$data_source, "geo")) {
  suppressPackageStartupMessages(library(GEOquery))

  if (is.null(cfg$geo$accession) || is.na(cfg$geo$accession)) {
    stop("data_source is 'geo' but config.yml's geo.accession is not set. ",
         "Fill it in once GEO assigns an accession after publication.")
  }

  dir.create(dirname(cfg$paths$raw_counts_csv), recursive = TRUE, showWarnings = FALSE)

  supp_files <- getGEOSuppFiles(cfg$geo$accession, baseDir = tempdir(), makeDirectory = TRUE)
  supp_path <- rownames(supp_files)[grepl(cfg$geo$supplementary_file, rownames(supp_files), fixed = TRUE)]
  stopifnot("Expected supplementary counts file not found in the GEO series" = length(supp_path) == 1)

  file.copy(supp_path, cfg$paths$raw_counts_csv, overwrite = TRUE)
  cat(sprintf("Fetched %s from GEO accession %s -> %s\n",
              cfg$geo$supplementary_file, cfg$geo$accession, cfg$paths$raw_counts_csv))
}
