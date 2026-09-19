#!/usr/bin/env Rscript
# scripts/09_build_supplement.R
# Assemble supplement/ (supplementary material for the ITEST 2026 paper,
# to be deposited on Zenodo) from stored pipeline results.
#
# Nothing is recomputed: the CSVs are copied from results/tables/ as-is.
# Rerun after scripts/run_all.R whenever results change.
#
# Usage (from the project root): Rscript scripts/09_build_supplement.R

proj_root <- here::here()
src_dir   <- file.path(proj_root, "results", "tables")
out_dir   <- file.path(proj_root, "supplement")
dir.create(out_dir, showWarnings = FALSE)

files <- c(
  "T1_theoretical.csv",
  "T2_monte_carlo.csv",
  "T6_misspecification.csv",
  "T7_bimodal_test.csv",
  "iris_versicolor_results.csv",
  "real_data_all_candidates.csv",
  "real_data_bootstrap.csv",
  "real_data_loo.csv"
)

missing <- files[!file.exists(file.path(src_dir, files))]
if (length(missing) > 0) stop("Missing in results/tables: ", paste(missing, collapse = ", "))

ok <- file.copy(file.path(src_dir, files), file.path(out_dir, files), overwrite = TRUE)
stopifnot(all(ok))
cat("Copied", length(files), "CSV files to", out_dir, "\n")

# Checksums let a reader confirm the deposit matches this repository.
md5 <- tools::md5sum(file.path(out_dir, files))
writeLines(paste(unname(md5), files), file.path(out_dir, "MD5SUMS.txt"))

# Package versions used to build the supplement.
pkgs <- c("EstemPMM", "ggplot2", "dplyr", "tidyr", "knitr", "patchwork")
ver  <- vapply(pkgs, function(p) {
  if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "not installed"
}, character(1))
writeLines(c(
  paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  R.version.string,
  paste("Platform:", R.version$platform),
  "",
  paste0(pkgs, " ", ver)
), file.path(out_dir, "SESSION_INFO.txt"))

cat("Wrote MD5SUMS.txt and SESSION_INFO.txt\n")
