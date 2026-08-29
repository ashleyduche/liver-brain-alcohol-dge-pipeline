suppressPackageStartupMessages({
  library(yaml)
})

cfg <- read_yaml("config.yml")

script_dir <- "R"
source(file.path(script_dir, "functions.R"))

cat(sprintf(
  "R %s | DESeq2 %s\n",
  paste(R.version$major, R.version$minor, sep = "."),
  as.character(packageVersion("DESeq2"))
))
