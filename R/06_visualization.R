# =============================================================================
# R/06_visualization.R
# ggplot2 visualization for all study blocks
# Requires: ggplot2, scales (loaded via .Rprofile / load_project)
# =============================================================================

#' Fig 1: Theoretical PMM3 efficiency curve g3(lambda) for TN distribution
plot_efficiency_curve <- function(T1_df) {
  library(ggplot2)
  ggplot(T1_df, aes(x = lambda, y = g3)) +
    geom_line(linewidth = 1.3, colour = "#2166AC") +
    geom_point(aes(colour = mode), size = 3.5) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    annotate("text", x = max(T1_df$lambda) * 0.6, y = 1.02,
             label = "OLS baseline (g₃ = 1)", colour = "grey40", size = 3.5) +
    annotate("vline", xintercept = 1, linetype = "dotted", colour = "orange") +
    scale_colour_manual(values = c(unimodal = "#4DAF4A", bimodal = "#E41A1C")) +
    scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                       limits = c(NA, 1.05)) +
    labs(
      title    = "Theoretical PMM3 Efficiency vs OLS for TN(λ)",
      subtitle = "g₃ < 1 means PMM3 has lower variance than OLS",
      x        = "Shape parameter λ",
      y        = "Variance ratio  g₃ = Var(PMM3) / Var(OLS)",
      colour   = "Distribution mode"
    ) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(colour = "grey40"))
}

#' Fig 2: Empirical vs theoretical g3 across sample sizes (Table T2 visual)
plot_g3_convergence <- function(mc_results, T1_df) {
  library(ggplot2)
  pmm3 <- mc_results[mc_results$method == "PMM3" & !is.na(mc_results$g3_empirical), ]
  pmm3 <- merge(pmm3, T1_df[, c("lambda", "g3")], by = "lambda", all.x = TRUE)

  ggplot(pmm3, aes(x = factor(n), y = g3_empirical)) +
    geom_col(fill = "#2166AC", alpha = 0.75, width = 0.6) +
    geom_hline(aes(yintercept = g3), colour = "#D73027",
               linetype = "dashed", linewidth = 1) +
    facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda))) +
    labs(
      title    = "Empirical g̃₃ vs Theoretical g₃ (PMM3 efficiency)",
      subtitle = "Red dashed = theoretical limit; bars = Monte Carlo estimate",
      x        = "Sample size n",
      y        = "g̃₃ = Var(PMM3) / Var(OLS)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.subtitle = element_text(colour = "grey40"))
}

#' Fig 3: ARE comparison across methods, lambda, and n
plot_are_comparison <- function(mc_results) {
  library(ggplot2)
  df <- mc_results[!is.na(mc_results$are), ]

  ggplot(df, aes(x = factor(n), y = are, fill = method)) +
    geom_col(position = "dodge", alpha = 0.85) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey30") +
    facet_wrap(~ lambda, labeller = label_bquote(lambda == .(lambda))) +
    scale_fill_brewer(palette = "Set1") +
    labs(
      title    = "Asymptotic Relative Efficiency vs OLS",
      subtitle = "ARE > 1: method more efficient than OLS",
      x        = "Sample size n",
      y        = "ARE = MSE(OLS) / MSE(method)",
      fill     = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(colour = "grey40"))
}

#' Fig 4: Bimodal test — convergence rate vs lambda
plot_bimodal_convergence <- function(bimodal_df) {
  library(ggplot2)
  ggplot(bimodal_df, aes(x = lambda, y = conv_rate,
                          colour = is_bimodal, shape = is_bimodal)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    geom_vline(xintercept = 1, linetype = "dotted", colour = "orange") +
    annotate("text", x = 1.05, y = 0.5, label = "bimodal boundary",
             colour = "orange", size = 3.5, hjust = 0) +
    scale_colour_manual(values = c("FALSE" = "#4DAF4A", "TRUE" = "#E41A1C"),
                        labels = c("Unimodal", "Bimodal")) +
    scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 17),
                       labels = c("Unimodal", "Bimodal")) +
    scale_y_continuous(limits = c(0, 1),
                       labels = function(x) paste0(x * 100, "%")) +
    labs(title  = "PMM3 Convergence Rate by λ (Bimodal Test)",
         x = "λ", y = "Convergence rate", colour = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "bottom")
}

#' Save a ggplot to results/figures/
save_plot <- function(p, filename, width = 10, height = 6.5, dpi = 300) {
  library(ggplot2)
  path <- file.path(CFG$paths$results_figures, filename)
  ggsave(path, plot = p, width = width, height = height, dpi = dpi)
  message("Saved: ", path)
  invisible(path)
}
