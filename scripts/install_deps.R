#!/usr/bin/env Rscript
# =============================================================================
# scripts/install_deps.R — Встановлення всіх залежностей проєкту
# Запускати один раз: Rscript scripts/install_deps.R
# =============================================================================

options(repos = c(CRAN = "https://cran.r-project.org"))

cat("=== PMM3-TN Research Project: Dependency Installation ===\n\n")

# CRAN пакети
cran_pkgs <- c(
  "ggplot2",    # Візуалізація
  "dplyr",      # Трансформації даних
  "tidyr",      # Reshape
  "knitr",      # Таблиці у Rmd
  "rmarkdown",  # Рендер HTML-звітів
  "nlme",       # Dataset nlme::Fatigue
  "pracma",     # Числові методи (hypergeo)
  "testthat",   # Unit-тести
  "scales",     # Форматування осей ggplot2
  "here",       # Управління робочою директорією
  "devtools"    # Для встановлення EstemPMM з GitHub
)

cat("--- Installing CRAN packages ---\n")
missing_cran <- cran_pkgs[!sapply(cran_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) {
  cat("Installing:", paste(missing_cran, collapse = ", "), "\n")
  install.packages(missing_cran, quiet = FALSE)
} else {
  cat("All CRAN packages already installed.\n")
}

# Перевірка після встановлення
cat("\n--- CRAN package status ---\n")
status <- sapply(cran_pkgs, requireNamespace, quietly = TRUE)
for (i in seq_along(status))
  cat(sprintf("  %-12s %s\n", names(status)[i], ifelse(status[i], "OK", "FAILED")))

# EstemPMM з GitHub
cat("\n--- Installing EstemPMM from GitHub ---\n")
if (!requireNamespace("EstemPMM", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    tryCatch({
      devtools::install_github("SZabolotnii/EstemPMM", quiet = FALSE)
      cat("EstemPMM installed:", requireNamespace("EstemPMM", quietly = TRUE), "\n")
    }, error = function(e) {
      cat("ERROR installing EstemPMM:", conditionMessage(e), "\n")
      cat("Try manually: devtools::install_github('SZabolotnii/EstemPMM')\n")
    })
  } else {
    cat("devtools not available. Install manually:\n")
    cat("  install.packages('devtools')\n")
    cat("  devtools::install_github('SZabolotnii/EstemPMM')\n")
  }
} else {
  cat("EstemPMM already installed.\n")
  # Перевірка наявності lm_pmm3
  ns       <- asNamespace("EstemPMM")
  has_pmm3 <- exists("lm_pmm3", envir = ns, inherits = FALSE)
  cat("lm_pmm3 in EstemPMM:", has_pmm3, "\n")
  if (!has_pmm3) {
    cat("NOTE: lm_pmm3 not found → pmm3_skeleton.R stub will be used.\n")
  }
}

cat("\n=== Installation complete ===\n")
cat("Next step: open RStudio, set working directory to project root,\n")
cat("then run: source('R/config.R'); load_project()\n")
