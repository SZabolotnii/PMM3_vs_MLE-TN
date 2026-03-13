# =============================================================================
# R/config.R — Централізована конфігурація дослідження
# PMM3 vs MLE-TN Comparative Analysis
# Salinas et al. (2023), Mathematics 11(5), 1271
# =============================================================================
# РЕДАГУЙ ЦЕЙ ФАЙЛ для зміни параметрів дослідження.
# Всі інші модулі читають CFG звідси.

CFG <- list(

  # --- Monte Carlo -----------------------------------------------------------
  # 500 для швидкого тесту, 5000 для фінального запуску (~30-60 хв)
  mc_replications = 500,
  sample_sizes    = c(46, 100, 200, 500),   # n=46 = biaxial fatigue
  lambda_grid     = c(0.5, 1.0, 1.5, 2.0, 3.0, 5.0),

  # --- Справжні параметри моделі (DGP) --------------------------------------
  beta_true = c(beta0 = 1.0, beta1 = 2.0),
  eta_true  = 1.0,

  # --- Порогові значення ----------------------------------------------------
  symmetry_threshold = 0.3,   # |gamma3| < threshold → симетричний
  gaussian_threshold = 0.1,   # |gamma4| < threshold → гаусовий

  # --- Розподіли для misspecification тесту (Блок 5) -----------------------
  misspec_distributions = c("uniform", "triangular", "logistic", "student10"),

  # --- Бімодальний тест lambda > 1 (Блок 6) --------------------------------
  bimodal_lambda_grid = c(0.8, 1.0, 1.2, 1.5, 2.0, 3.0),

  # --- Шляхи (відносні від PROJECT_ROOT) ------------------------------------
  paths = list(
    results_tables  = "results/tables",
    results_figures = "results/figures",
    mc_cache        = "results/mc_cache",
    data_raw        = "data/raw",
    data_processed  = "data/processed"
  ),

  # --- Відтворюваність ------------------------------------------------------
  seed         = 20260312,   # YYYYMMDD — дата початку дослідження
  use_parallel = FALSE,      # TRUE = parallel::mclapply (macOS safe)
  n_cores      = 2
)

# Функція завантаження всіх модулів проєкту
load_project <- function(verbose = TRUE) {
  r_files <- sort(list.files("R", pattern = "\\.R$", full.names = TRUE))
  r_files <- r_files[basename(r_files) != "config.R"]  # config вже завантажено
  for (f in r_files) {
    if (verbose) message("  Loading: ", basename(f))
    source(f, local = FALSE)
  }
  if (verbose) message("Project loaded. ", length(r_files), " modules.")
  invisible(r_files)
}

message(sprintf(
  "CFG loaded | M=%d | n=%s | lambda=%s",
  CFG$mc_replications,
  paste(CFG$sample_sizes, collapse = ","),
  paste(CFG$lambda_grid, collapse = ",")
))
