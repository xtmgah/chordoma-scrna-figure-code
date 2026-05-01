manifest <- read.delim("data_manifest.tsv", stringsAsFactors = FALSE)
manifest$exists <- file.exists(manifest$input_path)

missing_required <- manifest[manifest$required == "yes" & !manifest$exists, ]
missing_optional <- manifest[manifest$required != "yes" & !manifest$exists, ]

cat("Input file check\n")
cat("================\n")
cat("Required inputs:", sum(manifest$required == "yes"), "\n")
cat("Missing required:", nrow(missing_required), "\n")
cat("Missing optional:", nrow(missing_optional), "\n\n")

if (nrow(missing_required) > 0) {
  cat("Missing required files:\n")
  print(missing_required[, c("figure", "input_path", "notes")], row.names = FALSE)
  quit(status = 1)
}

cat("All required inputs are present.\n")
