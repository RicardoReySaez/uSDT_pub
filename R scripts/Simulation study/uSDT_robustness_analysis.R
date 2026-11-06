# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : uSDT_robustness_analysis.R                                 ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 29-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Assess hsdt robustness under Klauer (1998) and Goldstein (2022).
# H1: mean sensitivity difference; H2: true d' correlation.
# Goldstein defines its generating correlation on the accuracy scale.
# H3: hsdt intercept (beta_0) vs generating intercept (alpha), the mean
# influence sensitivity of unaware participants. These differ in Goldstein.
#
# Recompute bias, coverage and inference from the saved replications, as in
# uSDT_validation_analysis.R. Coverage uses available CIs; rejection uses
# testable fits, while rejection_all uses all replications. type_S is the
# significant wrong-sign rate.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(SimDesign)
library(dplyr)
library(tidyr)
library(ggplot2)
library(effectsize)
library(emmeans)
library(kableExtra)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

source("R scripts/R functions/plot_functions.R")

source("R scripts/R functions/simulation_functions.R")

dir.create("Results/Simulation studies/Processed", recursive = TRUE, showWarnings = FALSE)

N_colours <- c("30" = "#86b6ef", "100" = "#256abf", "200" = "#104281")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Prepare final data frame
# ─────────────────────────────────────────────────────────────────────────────

res  <- readRDS("Results/Simulation studies/Raw/Robustness study/uSDT_robustness_simres.rds")
reps <- SimResults(res)

# rho and intercept are generating values (intercept = alpha).
factors <- c("DGM", "N", "trials", "prop_unaware", "rho", "intercept")

# For H3, compare the same estimate against beta_0 and alpha.
est_col   <- c(diff = "diff", rho = "rho", intercept = "intercept",
               alpha = "intercept")
truth_col <- c(diff = "diff", rho = "cor_d", intercept = "intercept_hsdt",
               alpha = "intercept")

testable_only <- c("intercept", "alpha")

# Recovery and inference per design condition.
df_recovery <- bind_rows(lapply(names(est_col), function(p) {
  reps |>
    transmute(across(all_of(factors)),
              truth    = .data[[truth_col[[p]]]],
              est      = .data[[paste0(est_col[[p]], "_est")]],
              se       = .data[[paste0(est_col[[p]], "_se")]],
              lo       = .data[[paste0(est_col[[p]], "_lo")]],
              hi       = .data[[paste0(est_col[[p]], "_hi")]],
              pv       = .data[[paste0(est_col[[p]], "_p")]],
              testable = testable_H2 == 1,
              use      = failed == 0 & !is.na(est) &
                (testable | !p %in% testable_only)) |>
    summarise(truth         = first(truth),
              used          = mean(use),
              mean_est      = mean(est[use]),
              median_est    = median(est[use]),
              bias          = mean_est - truth,
              median_bias   = median_est - truth,
              RMSE          = sqrt(mean((est[use] - truth)^2)),
              emp_SD        = sd(est[use]),
              coverage      = mean((lo <= truth & hi >= truth)[use & !is.na(lo) & !is.na(hi)]),
              MSRSE         = mean(se[use & !is.na(se)]^2) / var(est[use & !is.na(se)]),
              missing_ci    = mean(is.na(lo[use]) | is.na(hi[use])),
              rejection     = mean(pv[!is.na(pv)] < .05),
              rejection_all = mean(!is.na(pv) & pv < .05),
              type_S        = if (abs(truth) > 1e-10) mean(!is.na(pv) & pv < .05 & sign(est) != sign(truth)) else NA_real_,
              .by = all_of(factors)) |>
    mutate(parameter = p)
}))

# Convergence and boundary measures from SimDesign.
df_boundary <- res |>
  transmute(across(all_of(factors)), failed, singular,
            untestable = untestable_H2, zero_sd_dD, zero_sd_dI, zero_sd_cD,
            one_rho, time_fit) |>
  mutate(parameter = "fit")

df_long <- bind_rows(
  pivot_longer(df_recovery, cols = -c(all_of(factors), parameter),
               names_to = "measure", values_to = "value"),
  pivot_longer(df_boundary, cols = -c(all_of(factors), parameter),
               names_to = "measure", values_to = "value"))

saveRDS(object = df_long, file = "Results/Simulation studies/Processed/robustness_results_longdata.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Convergence and boundary fits
# ─────────────────────────────────────────────────────────────────────────────

# Check convergence and boundary fits.
res |>
  summarise(failed = sum(failed * REPLICATIONS), converged = min(converged),
            .by = DGM)

all(res$untestable_H2 == res$untestable_H3)

df_long |>
  filter(parameter == "fit", measure != "time_fit") |>
  summarise(mean = mean(value), max = max(value), .by = c(DGM, measure)) |>
  as.data.frame()

df_long |>
  filter(parameter == "fit", measure == "untestable") |>
  summarise(value = mean(value), .by = c(DGM, N, trials)) |>
  pivot_wider(names_from = trials, names_prefix = "K_", values_from = value) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Global recovery
# ─────────────────────────────────────────────────────────────────────────────

# Goldstein: gap between the hsdt and generating intercepts.
res |>
  filter(DGM == "Goldstein", N == 30, trials == 10, intercept == 0) |>
  transmute(prop_unaware, rho, gap = intercept_hsdt - intercept) |>
  pivot_wider(names_from = rho, names_prefix = "rho_", values_from = gap) |>
  as.data.frame()

df_true_est <- df_long |>
  filter(parameter != "fit", measure %in% c("truth", "mean_est")) |>
  pivot_wider(names_from = measure, values_from = value)

# Concordance of condition means (CCC, precision and accuracy).
ccc_means <- df_true_est |>
  summarise(tibble::as_tibble_row(ccc(truth, mean_est)), .by = c(DGM, parameter))
ccc_means

# Concordance across individual replications.
ccc_reps <- bind_rows(lapply(names(est_col), function(p) {
  reps |>
    filter(failed == 0, !is.na(.data[[paste0(est_col[[p]], "_est")]]),
           testable_H2 == 1 | !p %in% testable_only) |>
    summarise(tibble::as_tibble_row(ccc(.data[[truth_col[[p]]]],
                                        .data[[paste0(est_col[[p]], "_est")]])),
              .by = c(DGM, N)) |>
    mutate(parameter = p)
}))
ccc_reps |>
  select(DGM, parameter, N, CCC) |>
  pivot_wider(names_from = N, names_prefix = "N_", values_from = CCC) |>
  as.data.frame()

# Bias that persists as N grows suggests model misspecification.
df_long |>
  filter(parameter != "fit", measure == "bias") |>
  summarise(value = mean(value), .by = c(DGM, parameter, N, trials)) |>
  pivot_wider(names_from = trials, names_prefix = "K_", values_from = value) |>
  arrange(DGM, parameter, N) |>
  as.data.frame()

# Check whether H1 bias comes from direct or indirect sensitivity.
res |>
  summarise(mu_dD = mean(mu_dD_bias), mu_dI = mean(mu_dI_bias),
            .by = c(DGM, trials)) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: ANOVA meta-models
# ─────────────────────────────────────────────────────────────────────────────

# Include up to three-way effects; higher orders are mostly Monte Carlo noise.
meta_formula <- value ~ (N + trials + prop_unaware + rho + intercept)^3

df_K <- filter(df_long, DGM == "Klauer")
df_G <- filter(df_long, DGM == "Goldstein")

# H2/H3 recovery and untestable fits (H3 evaluated against beta_0).
models_K <- list(
  B_H2    = aov(meta_formula, data = data_split(df_K, par = "rho", measure = "bias")),
  RMSE_H2 = aov(meta_formula, data = data_split(df_K, par = "rho", measure = "RMSE")),
  ECR_H2  = aov(meta_formula, data = data_split(df_K, par = "rho", measure = "coverage")),
  B_H3    = aov(meta_formula, data = data_split(df_K, par = "intercept", measure = "bias")),
  RMSE_H3 = aov(meta_formula, data = data_split(df_K, par = "intercept", measure = "RMSE")),
  ECR_H3  = aov(meta_formula, data = data_split(df_K, par = "intercept", measure = "coverage")),
  UNT     = aov(meta_formula, data = data_split(df_K, par = "fit", measure = "untestable")))

models_G <- list(
  B_H2    = aov(meta_formula, data = data_split(df_G, par = "rho", measure = "bias")),
  RMSE_H2 = aov(meta_formula, data = data_split(df_G, par = "rho", measure = "RMSE")),
  ECR_H2  = aov(meta_formula, data = data_split(df_G, par = "rho", measure = "coverage")),
  B_H3    = aov(meta_formula, data = data_split(df_G, par = "intercept", measure = "bias")),
  RMSE_H3 = aov(meta_formula, data = data_split(df_G, par = "intercept", measure = "RMSE")),
  ECR_H3  = aov(meta_formula, data = data_split(df_G, par = "intercept", measure = "coverage")),
  UNT     = aov(meta_formula, data = data_split(df_G, par = "fit", measure = "untestable")))

# Keep non-partial omega-squared effects >= .06 in any outcome.
ANOVA_es_K <- data.frame(
  Predictor = omega_squared(models_K$B_H2, partial = FALSE, ci = NULL)[, 1],
  sapply(models_K, function(fit)
    round(omega_squared(fit, partial = FALSE, ci = NULL)[, 2], 2))) |>
  filter(if_any(-Predictor, \(x) x >= .06))
ANOVA_es_G <- data.frame(
  Predictor = omega_squared(models_G$B_H2, partial = FALSE, ci = NULL)[, 1],
  sapply(models_G, function(fit)
    round(omega_squared(fit, partial = FALSE, ci = NULL)[, 2], 2))) |>
  filter(if_any(-Predictor, \(x) x >= .06))

labs_html  <- c(N = "N", trials = "K", prop_unaware = "π<sub>0</sub>",
                rho = "ρ", intercept = "α")
labs_latex <- c(N = "$N$", trials = "$K$", prop_unaware = "$\\pi_0$",
                rho = "$\\rho$", intercept = "$\\alpha$")

ANOVA_tbl_K <- ANOVA_es_K |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
                            paste(labs_html[x], collapse = " × ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("<strong>", txt_val, "</strong>"), txt_val)
         }))
ANOVA_tbl_G <- ANOVA_es_G |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
                            paste(labs_html[x], collapse = " × ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("<strong>", txt_val, "</strong>"), txt_val)
         }))

note_factors <- paste(
  "N = Sample size (30, 100 or 200). K = Trials per condition and task (10,",
  "30 or 50). π<sub>0</sub> = Proportion of unaware participants (.2, .5 or",
  ".8). ρ = Generating correlation (0 to .6). α = Generating intercept (0,",
  ".15 or .30). H3 is evaluated against the hsdt intercept, β<sub>0</sub>.",
  "RMSE = Root mean squared error.",
  "ECR = Empirical coverage rate of the 95% CI. UNT = Untestable H2 and H3.",
  "Effects with ω<sup>2</sup> < .06 in every measure are omitted, and",
  "effects equal or higher than .14 are marked in bold.")

# Klauer ANOVA table (HTML).
kbl(ANOVA_tbl_K,
    escape = FALSE,
    col.names = c("Predictor", rep(c("Bias", "RMSE", "ECR"), times = 2), "UNT"),
    align = c("l", rep("c", 7)),
    caption = "Effect sizes (non-partial ω<sup>2</sup>) with data generated from Klauer et al. (1998)") %>%
  add_header_above(c(" " = 1, "H2: Correlation" = 3, "H3: Intercept" = 3,
                     "Boundary" = 1), bold = TRUE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  footnote(general = note_factors, general_title = "Note.",
           footnote_as_chunk = TRUE, escape = FALSE)

# Goldstein ANOVA table (HTML).
kbl(ANOVA_tbl_G,
    escape = FALSE,
    col.names = c("Predictor", rep(c("Bias", "RMSE", "ECR"), times = 2), "UNT"),
    align = c("l", rep("c", 7)),
    caption = "Effect sizes (non-partial ω<sup>2</sup>) with data generated from Goldstein et al. (2022)") %>%
  add_header_above(c(" " = 1, "H2: Correlation" = 3, "H3: Intercept" = 3,
                     "Boundary" = 1), bold = TRUE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  footnote(general = note_factors, general_title = "Note.",
           footnote_as_chunk = TRUE, escape = FALSE)

# Export ANOVA tables for LaTeX.
ANOVA_tbl_K <- ANOVA_es_K |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
                            paste(labs_latex[x], collapse = " $\\times$ ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("\\textbf{", txt_val, "}"), txt_val)
         }))
ANOVA_tbl_G <- ANOVA_es_G |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
                            paste(labs_latex[x], collapse = " $\\times$ ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("\\textbf{", txt_val, "}"), txt_val)
         }))

save(ANOVA_tbl_K, ANOVA_tbl_G, ANOVA_es_K, ANOVA_es_G,
     file = "Results/Simulation studies/Processed/robustness_ANOVA_tables.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Estimated marginal means
# ─────────────────────────────────────────────────────────────────────────────

factor_html  <- c(N = "Sample size (N)", trials = "Trials per condition (K)",
                  prop_unaware = "Unaware participants (π<sub>0</sub>)",
                  rho = "Generating correlation (ρ)",
                  intercept = "Generating intercept (α)")
factor_latex <- c(N = "Sample size ($N$)", trials = "Trials per condition ($K$)",
                  prop_unaware = "Unaware participants ($\\pi_0$)",
                  rho = "Generating correlation ($\\rho$)",
                  intercept = "Generating intercept ($\\alpha$)")

# Bold coverages outside [.925, .975].
marginal_means_K <- cbind(
  marginal_means(anova_fit = models_K$B_H2, format = "html"),
  sapply(names(models_K)[-1], function(m)
    marginal_means(anova_fit = models_K[[m]], format = "html",
                   highlight_crit = if (grepl("ECR", m)) "coverage" else "none")$value))
marginal_means_G <- cbind(
  marginal_means(anova_fit = models_G$B_H2, format = "html"),
  sapply(names(models_G)[-1], function(m)
    marginal_means(anova_fit = models_G[[m]], format = "html",
                   highlight_crit = if (grepl("ECR", m)) "coverage" else "none")$value))

row_info <- table(marginal_means_K$Factor_Name)[unique(marginal_means_K$Factor_Name)]

# Marginal means tables (HTML).
marginal_means_K[, -1] %>%
  kable(
    format = "html",
    escape = FALSE,
    align = "lccccccc",
    col.names = c("", rep(c("Bias", "RMSE", "ECR"), times = 2), "UNT"),
    caption = "Estimated marginal means with data generated from Klauer et al. (1998)"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "H2: Correlation" = 3, "H3: Intercept" = 3,
                     "Boundary" = 1)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info)) {
      x <- pack_rows(x, factor_html[names(row_info)[i]], start_row,
                     start_row + row_info[i] - 1, escape = FALSE)
      start_row <- start_row + row_info[i]
    }
    x
  }}() %>%
  add_footnote(notation = "none", label = "Note: H3 is evaluated against the hsdt intercept, β<sub>0</sub>. Coverages outside [.925, .975] are in bold.")
marginal_means_G[, -1] %>%
  kable(
    format = "html",
    escape = FALSE,
    align = "lccccccc",
    col.names = c("", rep(c("Bias", "RMSE", "ECR"), times = 2), "UNT"),
    caption = "Estimated marginal means with data generated from Goldstein et al. (2022)"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "H2: Correlation" = 3, "H3: Intercept" = 3,
                     "Boundary" = 1)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info)) {
      x <- pack_rows(x, factor_html[names(row_info)[i]], start_row,
                     start_row + row_info[i] - 1, escape = FALSE)
      start_row <- start_row + row_info[i]
    }
    x
  }}() %>%
  add_footnote(notation = "none", label = "Note: H3 is evaluated against the hsdt intercept, β<sub>0</sub>. Coverages outside [.925, .975] are in bold.")

# Export marginal means for LaTeX.
marginal_means_K_latex <- cbind(
  marginal_means(anova_fit = models_K$B_H2, format = "latex"),
  sapply(names(models_K)[-1], function(m)
    marginal_means(anova_fit = models_K[[m]], format = "latex",
                   highlight_crit = if (grepl("ECR", m)) "coverage" else "none")$value))
marginal_means_G_latex <- cbind(
  marginal_means(anova_fit = models_G$B_H2, format = "latex"),
  sapply(names(models_G)[-1], function(m)
    marginal_means(anova_fit = models_G[[m]], format = "latex",
                   highlight_crit = if (grepl("ECR", m)) "coverage" else "none")$value))
marginal_means_K_latex$Factor_Name <- unname(factor_latex[marginal_means_K_latex$Factor_Name])
marginal_means_G_latex$Factor_Name <- unname(factor_latex[marginal_means_G_latex$Factor_Name])

save(marginal_means_K_latex, marginal_means_G_latex,
     file = "Results/Simulation studies/Processed/robustness_emmeans_tables.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Inference on the hypotheses
# ─────────────────────────────────────────────────────────────────────────────

# Summarise H2/H3 bias, coverage and null rejection by model.
# Flag liberal and conservative results using Bradley's limits.
scorecard <- df_long |>
  filter(parameter %in% c("rho", "intercept", "alpha"),
         measure %in% c("truth", "bias", "coverage", "rejection")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  mutate(null = abs(truth) < 1e-10) |>
  summarise(Bias     = sprintf("%.3f [%.3f, %.3f]", mean(bias), min(bias), max(bias)),
            ECR      = sprintf("%.3f [%.3f, %.3f]", mean(coverage), min(coverage), max(coverage)),
            ECR_lib  = sprintf("%.0f%%", 100 * mean(coverage < .925)),
            ECR_con  = sprintf("%.0f%%", 100 * mean(coverage > .975)),
            Null     = if (any(null)) sprintf("%.3f [%.3f]", mean(rejection[null]),
                                              max(rejection[null])) else "-",
            Null_lib = if (any(null)) sprintf("%.0f%%", 100 * mean(rejection[null] > .075)) else "-",
            Null_con = if (any(null)) sprintf("%.0f%%", 100 * mean(rejection[null] < .025)) else "-",
            .by = c(DGM, parameter)) |>
  arrange(factor(DGM, levels = c("Klauer", "Goldstein")),
          factor(parameter, levels = c("rho", "intercept", "alpha"))) |>
  mutate(parameter = recode(parameter, rho = "H2: Correlation",
                            intercept = "H3: Intercept (β<sub>0</sub>)", alpha = "H3: Intercept (α)"))

scorecard[, -1] %>%
  kable(
    format = "html",
    escape = FALSE,
    align = "lccccccc",
    col.names = c("", "Bias", "ECR", "< .925", "> .975", "Rejection under H0", "> .075", "< .025"),
    caption = "Robustness of the correlation and intercept hypotheses"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 3, "Liberal" = 1, "Conservative" = 1, " " = 1,
                     "Liberal" = 1, "Conservative" = 1)) %>%
  pack_rows("Klauer et al. (1998)", 1, 3) %>%
  pack_rows("Goldstein et al. (2022)", 4, 6) %>%
  add_footnote(notation = "none", label = "Note: Mean [minimum, maximum] over conditions for bias and ECR, and mean [maximum] over the null conditions for the rejection rate (testable fits). Liberal and conservative columns are percentages of conditions outside Bradley's intervals. α = generating intercept; β<sub>0</sub> = hsdt intercept of the population.")

# Export robustness summary for LaTeX.
scorecard_latex <- scorecard |>
  mutate(parameter = recode(parameter, "H3: Intercept (α)" = "H3: Intercept ($\\alpha$)",
                            "H3: Intercept (β<sub>0</sub>)" = "H3: Intercept ($\\beta_0$)"),
         across(everything(), \(x) gsub("%", "\\\\%", x)))
save(scorecard_latex, file = "Results/Simulation studies/Processed/robustness_summary_table.rdata")

# H3 false positives when alpha = 0.
df_long |>
  filter(parameter == "alpha", intercept == 0, measure == "rejection") |>
  summarise(value = mean(value), .by = c(DGM, prop_unaware, rho, N)) |>
  pivot_wider(names_from = N, names_prefix = "N_", values_from = value) |>
  arrange(DGM, prop_unaware, rho) |>
  as.data.frame()

# H2 false positives when rho = 0.
df_long |>
  filter(parameter == "rho", rho == 0, measure == "rejection") |>
  summarise(value = mean(value), .by = c(DGM, N, trials)) |>
  pivot_wider(names_from = trials, names_prefix = "K_", values_from = value) |>
  as.data.frame()

# H1 wrong-sign rate when the true difference is near zero.
df_long |>
  filter(parameter == "diff", measure %in% c("truth", "type_S")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  filter(abs(truth) < .1) |>
  summarise(mean = mean(type_S), max = max(type_S), conditions = n(), .by = DGM) |>
  as.data.frame()

# H2/H3 power among testable fits and all replications.
df_long |>
  filter((parameter == "rho" & rho > 0) | (parameter == "alpha" & intercept > 0),
         measure %in% c("rejection", "rejection_all")) |>
  summarise(value = mean(value), .by = c(DGM, parameter, measure, N)) |>
  pivot_wider(names_from = N, names_prefix = "N_", values_from = value) |>
  arrange(DGM, parameter, measure) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Figures of the supplemental material (Figures S9 to S11)
# ─────────────────────────────────────────────────────────────────────────────

dgm_levels <- c("Klauer", "Goldstein")
trials_axis <- scale_x_continuous(breaks = c(10, 30, 50))
panel_gap   <- theme(panel.spacing = unit(0.8, "lines"))

# Untestable H2/H3 fits.
fig_boundary <- df_long |>
  filter(parameter == "fit", measure == "untestable") |>
  summarise(value = mean(value), .by = c(DGM, N, trials, prop_unaware)) |>
  mutate(DGM    = factor(DGM, levels = dgm_levels),
         pi_lab = reorder(paste0("pi[0]=='", sub("^0\\.", ".", prop_unaware), "'"), prop_unaware)) |>
  ggplot(aes(x = trials, y = value, colour = factor(N))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_grid(DGM ~ pi_lab, labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  trials_axis +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Trials per condition", y = "Untestable H2 and H3") +
  theme_usdt() +
  panel_gap
save_figure(fig_boundary, "Supplemental/Figure_S9.png", width = 6.0, height = 3.6)

target_labels <- c(diff      = "mu[I]-mu[A]~(H1)",
                   rho       = "rho~(H2)",
                   intercept = "beta[0]~(H3)",
                   alpha     = "alpha~(H3)")

measure_labels <- c(bias = "Bias", coverage = "Coverage~of~the~'95%'~CI")
ref_lines <- data.frame(measure_lab = factor(measure_labels, levels = measure_labels),
                        y = c(0, .95))
ref_band  <- data.frame(measure_lab = factor(measure_labels["coverage"], levels = measure_labels),
                        ymin = .925, ymax = .975)

# Bias and coverage for rho, beta_0 and alpha.
fig_bias_coverage <- df_long |>
  filter(parameter %in% c("rho", "intercept", "alpha"), measure %in% c("bias", "coverage")) |>
  summarise(value = mean(value), .by = c(DGM, parameter, measure, N, trials)) |>
  mutate(DGM         = factor(DGM, levels = dgm_levels),
         target_lab  = factor(target_labels[parameter], levels = target_labels),
         measure_lab = factor(measure_labels[measure], levels = measure_labels)) |>
  ggplot(aes(x = trials, y = value, colour = factor(N))) +
  geom_rect(data = ref_band, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "grey92", inherit.aes = FALSE) +
  geom_hline(data = ref_lines, aes(yintercept = y), linetype = "22",
             colour = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_grid(measure_lab ~ DGM + target_lab, scales = "free_y", switch = "y",
             labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  trials_axis +
  labs(x = "Trials per condition", y = NULL) +
  theme_usdt() +
  theme(strip.placement = "outside") +
  panel_gap
save_figure(fig_bias_coverage, "Supplemental/Figure_S10.png", width = 6.9, height = 3.8)

# H3 false positives when alpha = 0.
fig_H3 <- df_long |>
  filter(parameter == "alpha", intercept == 0, measure == "rejection") |>
  summarise(value = mean(value), .by = c(DGM, N, rho, prop_unaware)) |>
  mutate(DGM    = factor(DGM, levels = dgm_levels),
         pi_lab = reorder(paste0("pi[0]=='", sub("^0\\.", ".", prop_unaware), "'"), prop_unaware)) |>
  ggplot(aes(x = rho, y = value, colour = factor(N))) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = .025, ymax = .075,
           fill = "grey92") +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey30",
             linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_grid(DGM ~ pi_lab, labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  scale_x_continuous(breaks = c(0, .2, .4, .6), labels = c("0", ".2", ".4", ".6")) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = expression("Generating correlation " * rho),
       y = expression("Rejection rate of H3 when " * alpha * " = 0")) +
  theme_usdt() +
  panel_gap
save_figure(fig_H3, "Supplemental/Figure_S11.png", width = 6.0, height = 3.6)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Reliability of the observed sensitivities
# ─────────────────────────────────────────────────────────────────────────────

# Reliability is the squared true-observed correlation in 100,000 draws.
# It depends on trials and, for awareness, the unaware proportion.
set.seed(2026)
reliab_rob <- expand.grid(DGM = c("Klauer", "Goldstein"), trials = c(10, 30, 50),
                          prop_unaware = c(.2, .5, .8), rho = c(0, .2, .4, .6),
                          intercept = c(0, .15, .30), stringsAsFactors = FALSE)
reliab_rob <- cbind(reliab_rob, t(mapply(robustness_reliability,
                                         DGM          = reliab_rob$DGM,
                                         prop_unaware = reliab_rob$prop_unaware,
                                         rho          = reliab_rob$rho,
                                         intercept    = reliab_rob$intercept,
                                         trials       = reliab_rob$trials)))

# Check that the validation design recovers nominal reliability.
reliab_check <- expand.grid(reliab = c(.1, .3, .5, .7), trials = c(10, 30, 50))
reliab_check <- cbind(reliab_check, t(mapply(function(reliab, trials) {
  pars <- population_parameters(reliab, trials, rho = .4, intercept = .15,
                                mu_dD = .2, mu_cD = 0, sd_cD = .3)
  d_D  <- rnorm(1e5, .2, pars[["sd_dD"]])
  d_I  <- pars[["mu_dI"]] + pars[["slope"]] * (d_D - .2) + rnorm(1e5, 0, pars[["sd_res"]])
  empirical_reliability(d_D, d_I, c_D = rnorm(1e5, 0, .3), k = trials)
}, reliab_check$reliab, reliab_check$trials)))
reliab_check
max(abs(c(reliab_check$reliab_D, reliab_check$reliab_I) - reliab_check$reliab))

reliab_rob |>
  summarise(direct = mean(reliab_D), indirect = mean(reliab_I),
            indirect_min = min(reliab_I), indirect_max = max(reliab_I),
            .by = c(DGM, prop_unaware, trials))

save(reliab_rob, reliab_check,
     file = "Results/Simulation studies/Processed/robustness_reliability.rdata")

# ─────────────────────────────────────────────────────────────────────────────
