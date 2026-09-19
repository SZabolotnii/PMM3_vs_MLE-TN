#!/usr/bin/env Rscript
# scripts/08_generate_paper_figures.R
# Generate all figures for the paper.
# Output: paper/figures/Fig{1-5}.png  (300 dpi, journal-ready)
#         paper/figures/Fig1_2_TN_densities_efficiency.png  (two-panel
#         figure used as Fig. 1 of the ITEST 2026 paper)
#
# Plots carry no in-plot title: Springer LNNS puts the caption under the
# figure, and long titles were clipped at the export width.
#
# Usage: Rscript scripts/08_generate_paper_figures.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
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

make_fig1 <- function(base_size = 12) {
ggplot(dens_df, aes(x = z, y = density, colour = lambda)) +
  geom_line(linewidth = 1.1) +
  geom_line(data = dnorm_df, aes(x = z, y = density),
            colour = "grey50", linewidth = 0.8, linetype = "dashed") +
  annotate("text", x = 0.3, y = dnorm(0) + 0.005, label = "N(0, 1)",
           colour = "grey40", size = base_size / 4, hjust = 0, vjust = 0) +
  # Skip the two lightest Blues shades: lambda = 0.5 was near-invisible.
  scale_colour_manual(values = RColorBrewer::brewer.pal(7, "Blues")[3:7]) +
  coord_cartesian(xlim = c(-4, 4), ylim = c(0, 0.45)) +
  labs(
    x       = "z",
    y       = "f(z | λ)",
    colour  = NULL
  ) +
  theme_bw(base_size = base_size) +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.88, 0.72),
        legend.key.height      = unit(0.7, "lines"),
        legend.key.width       = unit(1.2, "lines"),
        legend.margin     = margin(2, 4, 2, 4),
        legend.background = element_rect(fill = "white", colour = "grey80"))
}

save_fig(make_fig1(), "Fig1_TN_densities.png")

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

# TN(lambda) = equal mixture of N(-lambda, 1) and N(lambda, 1): bimodal iff
# lambda > 1 (see tn_cumulants()$is_bimodal). The earlier version drew the
# boundary at 1/sqrt(2), which contradicted T1 (lambda = 1 is unimodal).
bimodal_boundary <- 1

make_fig2 <- function(base_size = 12) {
ggplot(curve_df, aes(x = lambda, y = g3)) +
  geom_line(linewidth = 1.3, colour = "#2166AC") +
  geom_point(data = t1, aes(x = lambda, y = g3, colour = is_bimodal),
             size = 3.5, shape = 19) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = bimodal_boundary, linetype = "dotted",
             colour = "orange", linewidth = 0.9) +
  annotate("text", x = bimodal_boundary + 0.1, y = 0.15,
           label = "bimodal\nboundary", colour = "darkorange",
           size = base_size / 4, hjust = 0) +
  scale_colour_manual(values = c("FALSE" = "#4DAF4A", "TRUE" = "#E41A1C"),
                      labels = c("Unimodal", "Bimodal"),
                      name   = "TN mode") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2)) +
  labs(
    x       = "Shape parameter λ",
    y       = "g₃  (lower = more efficient)"
  ) +
  theme_bw(base_size = base_size) +
  theme(legend.position        = "inside",
        legend.position.inside = c(0.8, 0.75),
        legend.key.height      = unit(0.8, "lines"),
        legend.margin     = margin(2, 4, 2, 4),
        legend.background = element_rect(fill = "white", colour = "grey80"))
}

save_fig(make_fig2(), "Fig2_efficiency_curve.png")

# ── Fig 1+2: two-panel figure (paper Fig. 1) ─────────────────────────────────
cat("=== Fig 1+2: two-panel TN densities | g3 curve ===\n")

fig12 <- (make_fig1(base_size = 8) | make_fig2(base_size = 8)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

save_fig(fig12, "Fig1_2_TN_densities_efficiency.png", width = 24, height = 7)

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
  # One label per panel (the per-bar version overlapped the dashed line).
  geom_label(data = distinct(pmm3_mc, lam_label, g3),
             aes(x = 4.5, y = g3 + 0.08,
                 label = sprintf("theory %.3f", g3)),
             inherit.aes = FALSE, hjust = 1,
             colour = "#D73027", fill = "white", border.colour = NA,
             label.padding = unit(0.1, "lines"), size = 3.2) +
  facet_wrap(~ lam_label, nrow = 2) +
  labs(
    x        = "Sample size n",
    y        = "ĝ₃"
  ) +
  theme_bw(base_size = 11) +
  theme(strip.background = element_rect(fill = "#EEF4FB"))

# 18 x 9.8 cm keeps the aspect ratio the 10-page paper was laid out with
# (the old 18 x 11 export minus the clipped title band).
save_fig(fig3, "Fig3_g3_convergence.png", width = 18, height = 9.8)

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
    x        = "Sample size n",
    y        = "ARE",
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position   = "bottom",
        strip.background  = element_rect(fill = "#EEF4FB"))

# 18 x 10.8 cm: same reasoning as Fig 3.
save_fig(fig4, "Fig4_ARE_comparison.png", width = 18, height = 10.8)

# ── Fig 5: Iris versicolor — residual distribution ───────────────────────────
# Not used in the 10-page ITEST 2026 version (moved to the supplement).
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
      x        = "Residual",
      y        = "Density"
    ) +
    theme_bw(base_size = 12)
  save_fig(fig5, "Fig5_iris_residuals.png", width = 12, height = 9)
}

# Copy to paper/figures with standard name
cat("  Saved: Fig5_iris_residuals.png\n")

# ── Summary ───────────────────────────────────────────────────────────────────
cat("\n=== Done ===\n")
figs <- list.files(figures_out, "Fig.*\\.png")
cat("Generated", length(figs), "figures:\n")
cat(paste0("  ", figs, "\n"), sep = "")
