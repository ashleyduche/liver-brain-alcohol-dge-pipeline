suppressPackageStartupMessages({
  library(yaml)
})

# Pipeline settings
cfg <- read_yaml("config.yml")

# Shared helper functions and packages
script_dir <- "R"
source(file.path(script_dir, "functions.R"))

# Print the R and DESeq2 versions for the run log
cat(sprintf(
  "R %s | DESeq2 %s\n",
  paste(R.version$major, R.version$minor, sep = "."),
  as.character(packageVersion("DESeq2"))
))
