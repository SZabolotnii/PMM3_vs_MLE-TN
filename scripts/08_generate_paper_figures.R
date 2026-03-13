#!/usr/bin/env Rscript
# scripts/08_generate_paper_figures.R
# Generate all figures for the paper.
# Output: paper/figures/Fig{1-5}.png  (300 dpi, journal-ready)
#
# Usage: Rscript scripts/08_generate_paper_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

source("R/config.R")
load_project()

figures_out <- file.path(here::here(), "paper", "figures")
dir.create(figures_out, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(p, name, width = 14, height = 9, dpi = 300) {
  path <- file.path(figures_out, name)
  ggsave(path, plot = p, width = width, height = height,
         units = "cm", dpi = dpi)
  cat("  Saved:", name, "\n")
  invisible(path)
}

# ── Fig 1: TN density curves ──────────────────────────────────────────────────
cat("=== Fig 1: TN density curves ===\n")

lambdas <- c(0.5, 1.0, 1.5, 2.0, 3.0)
z_seq   <- seq(-4, 4, length.out = 400)

dens_df <- do.call(rbind, lapply(lambdas, function(lam) {
  d <- dtn(z_seq, xi = 0, eta = 1, lambda = lam)
  data.frame(z = z_seq, density = d,
             lambda = factor(paste0("λ = ", lam)))
}))

dnorm_df <- data.frame(z = z_seq, density = dnorm(z_seq),
                       lambda = factor("Normal (λ = 0)"))

fig1 <- ggplot(dens_df, aes(x = z, y = density, colour = lambda)) +
  geom_line(linewidth = 1.1) +
  geom_line(data = dnorm_df, aes(x = z, y = density),
            colour = "grey50", linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = 3.6, y = dnorm(0) * 0.06,
           label = "Normal", colour = "grey50", size = 3.2) +
  scale_colour_brewer(palette = "Blues", direction = 1) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 0.65)) +
  labs(
    title   = "Two-piece Normal (TN) density for selected shape parameters λ",
    x       = "z",
    y       = "f(z | λ)",
    colour  = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.82, 0.7),
        legend.background = element_rect(fill = "white", colour = "grey80"))

save_fig(fig1, "Fig1_TN_densities.png")

# ── Fig 2: Theoretical efficiency curve g3(lambda) ───────────────────────────
cat("=== Fig 2: Theoretical efficiency curve ===\n")

t1 <- read.csv("results/tables/T1_theoretical.csv")

# Dense lambda grid for smooth curve
lam_fine <- seq(0, 5, by = 0.05)
g3_fine  <- sapply(lam_fine, function(lam) {
  if (lam < 1e-6) return(1.0)
  cu     <- tn_cumulants(lam)
  denom  <- 6 + 9 * cu$gamma4 + cu$gamma6
  if (abs(denom) < 1e-10) return(1.0)
  1 - cu$gamma4^2 / denom
})
curve_df <- data.frame(lambda = lam_fine, g3 = g3_fine)

fig2 <- ggplot(curve_df, aes(x = lambda, y = g3)) +
  geom_line(linewidth = 1.3, colour = "#2166AC") +
  geom_point(data = t1, aes(x = lambda, y = g3, colour = is_bimodal),
             size = 3.5, shape = 19) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 1 / sqrt(2), linetype = "dotted",
             colour = "orange", linewidth = 0.9) +
  annotate("text", x = 1 / sqrt(2) + 0.1, y = 0.15,
           label = "bimodal\nboundary", colour = "darkorange",
           size = 3.0, hjust = 0) +
  scale_colour_manual(values = c("FALSE" = "#4DAF4A", "TRUE" = "#E41A1C"),
                      labels = c("Unimodal", "Bimodal"),
                      name   = "TN mode") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(
    title   = "Theoretical PMM3 efficiency  g₃ = Var(PMM3) / Var(OLS)",
    x       = "Shape parameter λ",
    y       = "g₃  (lower = more efficient)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill = "white", colour = "grey80"))

save_fig(fig2, "Fig2_efficiency_curve.png")

# ── Fig 3: Empirical g3 convergence to theory ────────────────────────────────
cat("=== Fig 3: Empirical g3 convergence ===\n")

mc <- read.csv("results/tables/T2_monte_carlo.csv")
t1_sub <- t1[, c("lambda", "g3")]

pmm3_mc <- mc |>
  filter(method == "PMM3", !is.na(g3_empirical)) |>
  merge(t1_sub, by = "lambda") |>
  mutate(
    lam_label = factor(paste0("λ = ", lambda),
                       levels = paste0("λ = ", sort(unique(lambda))))
  )

fig3 <- ggplot(pmm3_mc, aes(x = factor(n), y = g3_empirical)) +
  geom_col(fill = "#2166AC", alpha = 0.75, width = 0.6) +
  geom_hline(aes(yintercept = g3), colour = "#D73027",
             linetype = "dashed", linewidth = 1) +
  geom_text(aes(y = g3 + 0.04,
                label = sprintf("g₃=%.3f", g3)),
            colour = "#D73027", size = 2.8) +
  facet_wrap(~ lam_label, nrow = 2) +
  labs(
    title    = "Empirical g₃ vs theoretical (red dashed) — PMM3 efficiency",
    subtitle = "M = 500 replications; bars = ĝ₃ = Var(PMM3)/Var(OLS)",
    x        = "Sample size n",
    y        = "ĝ₃"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.subtitle = element_text(colour = "grey40"),
        strip.background = element_rect(fill = "#EEF4FB"))

save_fig(fig3, "Fig3_g3_convergence.png", width = 18, height = 11)

# ── Fig 4: ARE comparison (PMM3 vs MLE-TN vs OLS) ────────────────────────────
cat("=== Fig 4: ARE comparison ===\n")

are_df <- mc |>
  filter(method %in% c("MLE", "PMM3")) |>
  mutate(
    Method    = recode(method, "MLE" = "MLE-TN"),
    lam_label = factor(paste0("λ = ", lambda),
                       levels = paste0("λ = ", sort(unique(lambda))))
  )

fig4 <- ggplot(are_df, aes(x = factor(n), y = are, fill = Method)) +
  geom_col(position = position_dodge(0.7), width = 0.65, alpha = 0.88) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
  facet_wrap(~ lam_label, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("MLE-TN" = "#D73027", "PMM3" = "#2166AC")) +
  labs(
    title    = "ARE = Var(OLS)/Var(method) — Monte Carlo (M = 500)",
    subtitle = "ARE > 1: method more efficient than OLS; dashed = OLS baseline",
    x        = "Sample size n",
    y        = "ARE",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position   = "bottom",
        plot.subtitle     = element_text(colour = "grey40"),
        strip.background  = element_rect(fill = "#EEF4FB"))

save_fig(fig4, "Fig4_ARE_comparison.png", width = 18, height = 12)

# ── Fig 5: Iris versicolor — residual distribution ───────────────────────────
cat("=== Fig 5: Iris versicolor residual distribution ===\n")

# Load existing iris figure if available, else regenerate
iris_src <- "results/figures/iris_versicolor_resid_dist.png"
iris_dst <- file.path(figures_out, "Fig5_iris_residuals.png")

if (file.exists(iris_src)) {
  file.copy(iris_src, iris_dst, overwrite = TRUE)
  cat("  Copied from results/figures/\n")
} else {
  # Regenerate
  data(iris)
  vers <- subset(iris, Species == "versicolor")
  fit_ols <- lm(Sepal.Length ~ Sepal.Width, data = vers)
  resid_df <- data.frame(
    resid  = residuals(fit_ols),
    fitted = fitted(fit_ols)
  )
  fig5 <- ggplot(resid_df, aes(x = resid)) +
    geom_histogram(aes(y = after_stat(density)), bins = 12,
                   fill = "#2166AC", alpha = 0.7, colour = "white") +
    geom_density(linewidth = 1.1, colour = "#D73027") +
    stat_function(fun = dnorm,
                  args = list(mean = 0, sd = sd(resid_df$resid)),
                  linetype = "dashed", linewidth = 0.9, colour = "grey40") +
    labs(
      title    = "OLS residuals — Iris versicolor (n = 50)",
      subtitle = paste0("γ₄ = −0.966 (platykurtic) — PMM3 recommended"),
      x        = "Residual",
      y        = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))
  save_fig(fig5, "Fig5_iris_residuals.png", width = 12, height = 9)
}

# Copy to paper/figures with standard name
cat("  Saved: Fig5_iris_residuals.png\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Done ===\n")
figs <- list.files(figures_out, "Fig.*\\.png")
cat("Generated", length(figs), "figures:\n")
cat(paste0("  ", figs, "\n"), sep = "")
