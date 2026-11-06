# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : uSDT_validation_analysis.R                                 ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 29-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Analysis of the validation simulation study (uSDT_validation_sim.R):
# convergence and boundary fits, global recovery (true against mean estimated
# values, and Lin's concordance correlation), ANOVA meta-models of bias, RMSE
# and coverage, marginal means, and type I error and power of the three
# hypotheses.
#
# All the measures are recomputed from the stored replications, with these
# rules for boundary fits:
#   - Means, difference (H1), SDs and covariance: all the fits that did not
#     fail (a zero SD is a legitimate ML estimate).
#   - Correlation (H2): all the fits with an estimate, +-1 included.
#   - Intercept and slope (H3): testable fits only. In boundary fits the slope
#     divides the covariance by a zero variance, and reaches 1e38.
# Coverage and MSRSE use the fits with CI and SE. Rejection rates use the
# testable fits (rejection) or all the replications, an untestable hypothesis
# counting as not rejected (rejection_all).
#
# Effect sizes of the meta-models are non-partial omega-squared. Each cell is
# the mean of 1000 replications, so the residual is only Monte Carlo error and
# partial effect sizes approach 1 for most effects.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
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

# theme_usdt() and save_figure()
source("R scripts/R functions/plot_functions.R")

# data_split(), marginal_means() and ccc() (Section 4)
source("R scripts/R functions/simulation_functions.R")

# Output folder of the data and tables
dir.create("Results/Simulation studies/Processed", recursive = TRUE, showWarnings = FALSE)

# Colours of the sample sizes (one hue, light to dark)
N_colours <- c("30" = "#86b6ef", "60" = "#5598e7", "100" = "#256abf",
               "200" = "#104281")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Prepare final data frame
# ─────────────────────────────────────────────────────────────────────────────

# Load simulation results: summaries per condition, and every replication
res  <- readRDS("Results/Simulation studies/Raw/Validation study/uSDT_validation_simres.rds")
reps <- SimResults(res)

# Design factors and estimated parameters
factors <- c("N", "trials", "reliab", "rho", "intercept")
pars    <- c("mu_dD", "mu_dI", "mu_cD", "sd_dD", "sd_dI", "sd_cD", "cov_d",
             "diff", "rho", "intercept", "slope")

# Parameters only defined in fits with testable hypotheses
testable_only <- c("intercept", "slope")

# Recovery of each parameter per condition
df_recovery <- bind_rows(lapply(pars, function(p) {
  reps |>
    transmute(across(all_of(factors)),
              truth    = .data[[p]],
              est      = .data[[paste0(p, "_est")]],
              se       = .data[[paste0(p, "_se")]],
              lo       = .data[[paste0(p, "_lo")]],
              hi       = .data[[paste0(p, "_hi")]],
              pv       = .data[[paste0(p, "_p")]],
              testable = testable_H2 == 1,
              use      = failed == 0 & !is.na(est) &
                         (testable | !p %in% testable_only)) |>
    summarise(truth         = first(truth),
              used          = mean(use),
              mean_est      = mean(est[use]),
              median_est    = median(est[use]),
              bias          = mean_est - truth,
              rel_bias      = if (truth != 0) bias / truth else NA_real_,
              median_bias   = median_est - truth,
              bias_testable = mean(est[use & testable]) - truth,
              RMSE          = sqrt(mean((est[use] - truth)^2)),
              emp_SD        = sd(est[use]),
              coverage      = mean((lo <= truth & hi >= truth)[use & !is.na(lo) & !is.na(hi)]),
              MSRSE         = mean(se[use & !is.na(se)]^2) / var(est[use & !is.na(se)]),
              missing_ci    = mean(is.na(lo[use]) | is.na(hi[use])),
              rejection     = if (any(!is.na(pv))) mean(pv[!is.na(pv)] < .05) else NA_real_,
              rejection_all = if (any(!is.na(pv))) mean(!is.na(pv) & pv < .05) else NA_real_,
              .by = all_of(factors)) |>
    mutate(parameter = p)
}))

# Fits and boundary, from the summaries of SimDesign: failed and singular fits,
# untestable H2 and H3 (always together), zero SDs and correlation at +-1
df_boundary <- res |>
  transmute(across(all_of(factors)), failed, singular,
            untestable = untestable_H2, zero_sd_dD, zero_sd_dI, zero_sd_cD,
            one_rho, time_fit) |>
  mutate(parameter = "fit")

# Long data frame: condition, parameter, measure and value
df_long <- bind_rows(
  pivot_longer(df_recovery, cols = -c(all_of(factors), parameter),
               names_to = "measure", values_to = "value"),
  pivot_longer(df_boundary, cols = -c(all_of(factors), parameter),
               names_to = "measure", values_to = "value"))

# Save this final df
saveRDS(object = df_long, file = "Results/Simulation studies/Processed/validation_results_longdata.rds")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Convergence and boundary fits
# ─────────────────────────────────────────────────────────────────────────────

# Failed fits (hsdt() returned an error) and convergence of the rest
sum(res$failed * res$REPLICATIONS)
range(res$converged)

# H1 is always testable; H2 and H3 are testable or untestable together
max(res$untestable_H1)
all(res$untestable_H2 == res$untestable_H3)

# Every singular fit is untestable, and a few non-singular fits too
table(singular = reps$singular, testable = reps$testable_H2, useNA = "ifany")

# Mean and maximum of each boundary measure over conditions
df_long |>
  filter(parameter == "fit", measure != "time_fit") |>
  summarise(mean = mean(value), max = max(value), .by = measure) |>
  as.data.frame()

# Untestable fits by sample size and reliability (mean over the other factors)
df_long |>
  filter(parameter == "fit", measure == "untestable") |>
  summarise(value = mean(value), .by = c(N, reliab)) |>
  pivot_wider(names_from = reliab, values_from = value) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Global recovery
# ─────────────────────────────────────────────────────────────────────────────

# Parameters whose true value varies across conditions (mu_dD, mu_cD and sd_cD
# are fixed, so they are summarised by their bias, RMSE and coverage)
pars_vary <- c("mu_dI", "diff", "sd_dD", "sd_dI", "cov_d", "rho", "intercept",
               "slope")

# Parameters with a fixed true value
df_long |>
  filter(parameter %in% c("mu_dD", "mu_cD", "sd_cD"),
         measure %in% c("bias", "RMSE", "coverage", "MSRSE")) |>
  summarise(min = min(value), max = max(value), .by = c(parameter, measure)) |>
  as.data.frame()

# Mean estimate and true value per condition
df_true_est <- df_long |>
  filter(parameter %in% pars_vary, measure %in% c("truth", "mean_est")) |>
  pivot_wider(names_from = measure, values_from = value)

# Concordance of the mean estimates with the truth: CCC, and its precision (r)
# and accuracy (C_b) parts
ccc_means <- df_true_est |>
  summarise(tibble::as_tibble_row(ccc(truth, mean_est)), .by = parameter)
ccc_means

# Concordance of the mean estimates by reliability. Only comparable across
# reliabilities for the means and hypotheses: the true SDs and covariance
# change with reliability, and so does their range
ccc_reliab <- df_true_est |>
  summarise(tibble::as_tibble_row(ccc(truth, mean_est)), .by = c(parameter, reliab))
ccc_reliab |>
  select(parameter, reliab, CCC) |>
  pivot_wider(names_from = reliab, values_from = CCC) |>
  as.data.frame()

# Concordance of single-replication estimates with the truth, by sample size
ccc_reps <- bind_rows(lapply(pars_vary, function(p) {
  reps |>
    filter(failed == 0, !is.na(.data[[paste0(p, "_est")]]),
           testable_H2 == 1 | !p %in% testable_only) |>
    summarise(tibble::as_tibble_row(ccc(.data[[p]], .data[[paste0(p, "_est")]])),
              .by = N) |>
    mutate(parameter = p)
}))
ccc_reps |>
  select(parameter, N, CCC) |>
  pivot_wider(names_from = N, names_prefix = "N_", values_from = CCC) |>
  as.data.frame()

# Intercept and slope have a few extreme estimates even in testable fits: the
# median bias is a robust check of the mean bias
df_long |>
  filter(parameter %in% testable_only, measure %in% c("bias", "median_bias")) |>
  summarise(min = min(value), median = median(value), max = max(value),
            .by = c(parameter, measure)) |>
  as.data.frame()

# Correlation: bias over all fits and over testable fits (selection effect)
df_long |>
  filter(parameter == "rho", measure %in% c("bias", "bias_testable")) |>
  summarise(value = mean(value), .by = c(measure, reliab, rho)) |>
  pivot_wider(names_from = rho, names_prefix = "rho_", values_from = value) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: ANOVA meta-models
# ─────────────────────────────────────────────────────────────────────────────

# Awareness-task SD meta-models: relative bias, RMSE and coverage
RB_sdD   <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dD", measure = "rel_bias"))
RMSE_sdD <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dD", measure = "RMSE"))
ECR_sdD  <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dD", measure = "coverage"))

# Influence-task SD meta-models: relative bias, RMSE and coverage
RB_sdI   <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dI", measure = "rel_bias"))
RMSE_sdI <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dI", measure = "RMSE"))
ECR_sdI  <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "sd_dI", measure = "coverage"))

# Boundary meta-model: untestable H2 and H3
UNT <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
           data = data_split(df_long = df_long, par = "fit", measure = "untestable"))

# H1 meta-models: bias, RMSE and coverage of the difference
B_diff    <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                 data = data_split(df_long = df_long, par = "diff", measure = "bias"))
RMSE_diff <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                 data = data_split(df_long = df_long, par = "diff", measure = "RMSE"))
ECR_diff  <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                 data = data_split(df_long = df_long, par = "diff", measure = "coverage"))

# H2 meta-models: bias, RMSE and coverage of the correlation
B_rho    <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "rho", measure = "bias"))
RMSE_rho <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "rho", measure = "RMSE"))
ECR_rho  <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "rho", measure = "coverage"))

# H3 meta-models: bias, RMSE and coverage of the intercept
B_int    <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "intercept", measure = "bias"))
RMSE_int <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "intercept", measure = "RMSE"))
ECR_int  <- aov(value ~ (N + trials + reliab + rho + intercept)^4,
                data = data_split(df_long = df_long, par = "intercept", measure = "coverage"))

# Meta-models of each table
models_re  <- list(RB_D = RB_sdD, RMSE_D = RMSE_sdD, ECR_D = ECR_sdD,
                   RB_I = RB_sdI, RMSE_I = RMSE_sdI, ECR_I = ECR_sdI, UNT = UNT)
models_hyp <- list(B_H1 = B_diff, RMSE_H1 = RMSE_diff, ECR_H1 = ECR_diff,
                   B_H2 = B_rho,  RMSE_H2 = RMSE_rho,  ECR_H2 = ECR_rho,
                   B_H3 = B_int,  RMSE_H3 = RMSE_int,  ECR_H3 = ECR_int)

# ANOVA results based on effect sizes (non-partial omega-squared), keeping the
# effects of at least .06 (medium) in one measure: smaller high-order effects
# are mostly Monte Carlo error
ANOVA_es_re <- data.frame(
  Predictor = omega_squared(RB_sdD, partial = FALSE, ci = NULL)[, 1],
  sapply(models_re, function(fit)
    round(omega_squared(fit, partial = FALSE, ci = NULL)[, 2], 2))) |>
  filter(if_any(-Predictor, \(x) x >= .06))
ANOVA_es_hyp <- data.frame(
  Predictor = omega_squared(B_diff, partial = FALSE, ci = NULL)[, 1],
  sapply(models_hyp, function(fit)
    round(omega_squared(fit, partial = FALSE, ci = NULL)[, 2], 2))) |>
  filter(if_any(-Predictor, \(x) x >= .06))

# Predictor labels (K = trials per condition, R = reliability)
labs_html  <- c(N = "N", trials = "K", reliab = "R", rho = "ρ",
                intercept = "β<sub>0</sub>")
labs_latex <- c(N = "$N$", trials = "$K$", reliab = "$R$", rho = "$\\rho$",
                intercept = "$\\beta_0$")

# Set in bold omega-squared values above or equal to 0.14
ANOVA_tbl_re <- ANOVA_es_re |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
           paste(labs_html[x], collapse = " × ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("<strong>", txt_val, "</strong>"), txt_val)
         }))
ANOVA_tbl_hyp <- ANOVA_es_hyp |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
           paste(labs_html[x], collapse = " × ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("<strong>", txt_val, "</strong>"), txt_val)
         }))

# Common note of the ANOVA tables
note_factors <- paste(
  "N = Sample size (30, 60, 100 or 200). K = Trials per condition and task",
  "(10, 30 or 50). R = Reliability (.1 to .7). ρ = True correlation",
  "(0 to .8). β<sub>0</sub> = True intercept (0, .15 or .30). Effects with",
  "ω<sup>2</sup> < .06 in every measure are omitted, and effects equal",
  "or higher than .14 are marked in bold.")

# Final table results: random-effect SDs and boundary
kbl(ANOVA_tbl_re,
    escape = FALSE,
    col.names = c("Predictor", "RB", "RMSE", "ECR", "RB", "RMSE", "ECR", "UNT"),
    align = c("l", rep("c", 7)),
    caption = "Effect sizes (non-partial ω<sup>2</sup>) for the random-effect SDs and boundary fits") %>%
  add_header_above(c(" " = 1, "Awareness SD" = 3, "Influence SD" = 3,
                     "Boundary" = 1), bold = TRUE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  footnote(general = paste(note_factors, "RB = Relative bias. RMSE = Root mean squared error. ECR = Empirical coverage rate of the 95% CI. UNT = Untestable H2 and H3."),
           general_title = "Note.", footnote_as_chunk = TRUE, escape = FALSE)

# Final table results: hypotheses
kbl(ANOVA_tbl_hyp,
    escape = FALSE,
    col.names = c("Predictor", rep(c("Bias", "RMSE", "ECR"), times = 3)),
    align = c("l", rep("c", 9)),
    caption = "Effect sizes (non-partial ω<sup>2</sup>) for the three hypotheses") %>%
  add_header_above(c(" " = 1, "H1: Difference" = 3, "H2: Correlation" = 3,
                     "H3: Intercept" = 3), bold = TRUE) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  footnote(general = paste(note_factors, "RMSE = Root mean squared error. ECR = Empirical coverage rate of the 95% CI."),
           general_title = "Note.", footnote_as_chunk = TRUE, escape = FALSE)

# Prepare the same ANOVA tables but just for LaTeX and save them
ANOVA_tbl_re <- ANOVA_es_re |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
           paste(labs_latex[x], collapse = " $\\times$ ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("\\textbf{", txt_val, "}"), txt_val)
         }))
ANOVA_tbl_hyp <- ANOVA_es_hyp |>
  mutate(Predictor = sapply(strsplit(Predictor, ":"), \(x)
           paste(labs_latex[x], collapse = " $\\times$ ")),
         across(where(is.numeric), function(x) {
           txt_val <- sprintf("%.2f", x)
           ifelse(x >= 0.14, paste0("\\textbf{", txt_val, "}"), txt_val)
         }))

save(ANOVA_tbl_re, ANOVA_tbl_hyp, ANOVA_es_re, ANOVA_es_hyp,
     file = "Results/Simulation studies/Processed/validation_ANOVA_tables.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Estimated marginal means
# ─────────────────────────────────────────────────────────────────────────────

# Labels of the design factors
factor_html  <- c(N = "Sample size (N)", trials = "Trials per condition (K)",
                  reliab = "Reliability (R)", rho = "True correlation (ρ)",
                  intercept = "True intercept (β<sub>0</sub>)")
factor_latex <- c(N = "Sample size ($N$)", trials = "Trials per condition ($K$)",
                  reliab = "Reliability ($R$)", rho = "True correlation ($\\rho$)",
                  intercept = "True intercept ($\\beta_0$)")

# Random-effect SDs and boundary
marginal_means_re <- cbind(
  marginal_means(anova_fit = RB_sdD, highlight_crit = "relative", format = "html"),
  marginal_means(anova_fit = RMSE_sdD, format = "html")[-c(1:2)],
  marginal_means(anova_fit = ECR_sdD, highlight_crit = "coverage", format = "html")[-c(1:2)],
  marginal_means(anova_fit = RB_sdI, highlight_crit = "relative", format = "html")[-c(1:2)],
  marginal_means(anova_fit = RMSE_sdI, format = "html")[-c(1:2)],
  marginal_means(anova_fit = ECR_sdI, highlight_crit = "coverage", format = "html")[-c(1:2)],
  marginal_means(anova_fit = UNT, format = "html")[-c(1:2)])

# Compute row grouping (unique preserves factor appearance order)
row_info <- table(marginal_means_re$Factor_Name)[unique(marginal_means_re$Factor_Name)]

# Generate HTML table (select by position to avoid duplicate column name errors)
marginal_means_re[, -1] %>%
  kable(
    format = "html",
    escape = FALSE,
    align = "lccccccc",
    col.names = c("", "RB", "RMSE", "ECR", "RB", "RMSE", "ECR", "UNT"),
    caption = "Estimated marginal means for the random-effect SDs and boundary fits"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "Awareness SD" = 3, "Influence SD" = 3, "Boundary" = 1)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info)) {
      x <- pack_rows(x, factor_html[names(row_info)[i]], start_row,
                     start_row + row_info[i] - 1, escape = FALSE)
      start_row <- start_row + row_info[i]
    }
    x
  }}() %>%
  add_footnote(notation = "none", label = "Note: Relative biases of 10% or more, and coverages outside [.925, .975], are in bold.")

# Hypotheses
marginal_means_hyp <- cbind(
  marginal_means(anova_fit = B_diff, format = "html"),
  marginal_means(anova_fit = RMSE_diff, format = "html")[-c(1:2)],
  marginal_means(anova_fit = ECR_diff, highlight_crit = "coverage", format = "html")[-c(1:2)],
  marginal_means(anova_fit = B_rho, format = "html")[-c(1:2)],
  marginal_means(anova_fit = RMSE_rho, format = "html")[-c(1:2)],
  marginal_means(anova_fit = ECR_rho, highlight_crit = "coverage", format = "html")[-c(1:2)],
  marginal_means(anova_fit = B_int, format = "html")[-c(1:2)],
  marginal_means(anova_fit = RMSE_int, format = "html")[-c(1:2)],
  marginal_means(anova_fit = ECR_int, highlight_crit = "coverage", format = "html")[-c(1:2)])

# Generate HTML table
marginal_means_hyp[, -1] %>%
  kable(
    format = "html",
    escape = FALSE,
    align = "lccccccccc",
    col.names = c("", rep(c("Bias", "RMSE", "ECR"), times = 3)),
    caption = "Estimated marginal means for the three hypotheses"
  ) %>%
  kable_classic_2(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "H1: Difference" = 3, "H2: Correlation" = 3,
                     "H3: Intercept" = 3)) %>%
  { \(x) {
    start_row <- 1
    for(i in seq_along(row_info)) {
      x <- pack_rows(x, factor_html[names(row_info)[i]], start_row,
                     start_row + row_info[i] - 1, escape = FALSE)
      start_row <- start_row + row_info[i]
    }
    x
  }}() %>%
  add_footnote(notation = "none", label = "Note: Coverages outside [.925, .975] are in bold.")

# Save estimated marginal means in LaTeX format
marginal_means_re_latex <- cbind(
  marginal_means(anova_fit = RB_sdD, highlight_crit = "relative", format = "latex"),
  marginal_means(anova_fit = RMSE_sdD, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = ECR_sdD, highlight_crit = "coverage", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = RB_sdI, highlight_crit = "relative", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = RMSE_sdI, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = ECR_sdI, highlight_crit = "coverage", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = UNT, format = "latex")[-c(1:2)])
marginal_means_hyp_latex <- cbind(
  marginal_means(anova_fit = B_diff, format = "latex"),
  marginal_means(anova_fit = RMSE_diff, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = ECR_diff, highlight_crit = "coverage", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = B_rho, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = RMSE_rho, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = ECR_rho, highlight_crit = "coverage", format = "latex")[-c(1:2)],
  marginal_means(anova_fit = B_int, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = RMSE_int, format = "latex")[-c(1:2)],
  marginal_means(anova_fit = ECR_int, highlight_crit = "coverage", format = "latex")[-c(1:2)])

# Factor labels of the LaTeX tables (duplicated column names, so base R)
marginal_means_re_latex$Factor_Name  <- unname(factor_latex[marginal_means_re_latex$Factor_Name])
marginal_means_hyp_latex$Factor_Name <- unname(factor_latex[marginal_means_hyp_latex$Factor_Name])

save(marginal_means_re_latex, marginal_means_hyp_latex,
     file = "Results/Simulation studies/Processed/validation_emmeans_tables.rdata")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Type I error and power of the hypotheses
# ─────────────────────────────────────────────────────────────────────────────

# H1 has no null condition: the smallest true difference is
min(abs(res$diff))

# Type I error (H2: rho = 0; H3: intercept = 0, and slope when rho = 0), and
# proportion of conditions inside Bradley's liberal interval [.025, .075]
bind_rows(
  df_long |> filter(parameter == "rho", rho == 0),
  df_long |> filter(parameter == "intercept", intercept == 0),
  df_long |> filter(parameter == "slope", rho == 0)) |>
  filter(measure %in% c("rejection", "rejection_all")) |>
  summarise(mean = mean(value), min = min(value), max = max(value),
            bradley = mean(value >= .025 & value <= .075),
            .by = c(parameter, measure)) |>
  as.data.frame()

# Type I error of H2 over the testable fits, by sample size and reliability
df_long |>
  filter(parameter == "rho", rho == 0, measure == "rejection") |>
  summarise(value = mean(value), .by = c(N, reliab)) |>
  pivot_wider(names_from = reliab, values_from = value) |>
  as.data.frame()

# Conditions where the type I error of H3 (intercept) is above .075
df_long |>
  filter(parameter == "intercept", intercept == 0, measure == "rejection",
         value > .075) |>
  select(all_of(factors), value) |>
  as.data.frame()

# Power of H2 and H3 over the testable fits and over all replications, by
# sample size and reliability (mean over the other factors)
df_long |>
  filter((parameter == "rho" & rho > 0) | (parameter == "intercept" & intercept > 0),
         measure %in% c("rejection", "rejection_all")) |>
  summarise(value = mean(value), .by = c(parameter, measure, N, reliab)) |>
  pivot_wider(names_from = reliab, values_from = value) |>
  arrange(parameter, measure, N) |>
  as.data.frame()

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Figures of the manuscript (Figures 6 and 7)
# ─────────────────────────────────────────────────────────────────────────────

# Labels of the population parameters and the intercept, in the order of the
# columns of the figure: means, SDs, and correlation and intercept. The
# awareness-task mean, fixed at .2, is replaced by the group-level difference
ms_labels <- c(mu_dI = "Influence~mean~mu[I]", diff = "Difference~mu[I]-mu[A]",
               sd_dI = "Influence~SD~sigma[I]", sd_dD = "Awareness~SD~sigma[A]",
               rho = "Correlation~rho", intercept = "Intercept~beta[0]")

# True value against mean estimate per condition, with the CCC of each
# parameter
ccc_ms_labels <- ccc_means |>
  filter(parameter %in% names(ms_labels)) |>
  mutate(panel = factor(ms_labels[parameter], levels = ms_labels),
         label = paste("CCC =", sub("^0\\.", ".", sprintf("%.3f", CCC))))
fig_sim1_recovery <- df_long |>
  filter(parameter %in% names(ms_labels), measure %in% c("truth", "mean_est")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  mutate(panel = factor(ms_labels[parameter], levels = ms_labels)) |>
  arrange(reliab) |>
  ggplot(aes(x = truth, y = mean_est, colour = reliab)) +
  geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey30",
              linewidth = 0.5) +
  geom_blank(aes(x = mean_est, y = truth)) +
  geom_point(size = 0.7) +
  geom_text(data = ccc_ms_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.08, vjust = 1.5, size = 2.4, inherit.aes = FALSE) +
  facet_wrap(~ panel, nrow = 2, dir = "v", scales = "free", labeller = label_parsed) +
  scale_colour_gradient(low = "#86b6ef", high = "#0d366b", breaks = c(.1, .4, .7),
                        labels = c(".1", ".4", ".7"), name = "Reliability") +
  labs(x = "True value", y = "Mean estimate") +
  theme_usdt() +
  theme(legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.25, "cm"),
        plot.margin = margin(5.5, 5.5, 14, 5.5))
save_figure(fig_sim1_recovery, "Manuscript/Figure_6.png", width = 6.6, height = 4.0)

# Power of H2 and H3 by reliability, sample size and effect size, averaged
# over trials and the other parameter, over all replications (untestable fits
# count as not rejected; manuscript) or over the testable fits (supplemental
# material). Cells with power of .80 or more in bold
fig_power <- lapply(c(all = "rejection_all", testable = "rejection"), function(m) {
  df_power <- df_long |>
    filter(measure == m,
           (parameter == "rho" & rho > 0) | (parameter == "intercept" & intercept > 0)) |>
    mutate(effect = ifelse(parameter == "rho", rho, intercept)) |>
    summarise(power = mean(value), .by = c(parameter, effect, N, reliab)) |>
    arrange(parameter != "rho", effect) |>
    mutate(panel = ifelse(parameter == "rho",
                          paste0("H2:~rho=='", sub("^0\\.", ".", effect), "'"),
                          paste0("H3:~beta[0]=='", sub("^0\\.", ".", sprintf("%.2f", effect)), "'")),
           panel = factor(panel, levels = unique(panel)),
           label = ifelse(power > .995, "1", sub("^0\\.", ".", sprintf("%.2f", power))))
  ggplot(df_power, aes(x = factor(reliab), y = factor(N), fill = power)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_text(aes(label = label, colour = power > .55,
                  fontface = ifelse(power >= .80, "bold", "plain")), size = 2.3) +
    facet_wrap(~ panel, ncol = 3, labeller = label_parsed) +
    scale_x_discrete(labels = \(x) sub("^0\\.", ".", x)) +
    scale_fill_gradientn(colours = c("#eef4fc", "#86b6ef", "#256abf", "#0d366b"),
                         limits = c(0, 1), breaks = c(0, .5, .8, 1),
                         labels = c("0", ".50", ".80", "1"), name = "Power") +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
    labs(x = "Reliability (R)", y = "Sample size (N)") +
    theme_usdt() +
    theme(panel.border = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(),
          legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.25, "cm"),
          plot.margin = margin(5.5, 5.5, 14, 5.5))
})
save_figure(fig_power$all, "Manuscript/Figure_7.png", width = 6.9, height = 4.4)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 11: Figures of the supplemental material (Figures S1 to S8)
# ─────────────────────────────────────────────────────────────────────────────

# Reliability axis, and facet labels of the true values, without leading zeros
reliab_axis <- scale_x_continuous(breaks = c(.1, .4, .7),
                                  labels = c(".1", ".4", ".7"))
panel_gap   <- theme(panel.spacing = unit(0.8, "lines"))

# Untestable H2 and H3 by reliability and sample size, per true correlation
fig_boundary <- df_long |>
  filter(parameter == "fit", measure == "untestable") |>
  summarise(value = mean(value), .by = c(N, reliab, rho)) |>
  mutate(rho_lab = reorder(paste0("rho=='", sub("^0\\.", ".", rho), "'"), rho)) |>
  ggplot(aes(x = reliab, y = value, colour = factor(N))) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_wrap(~ rho_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  reliab_axis +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Reliability", y = "Untestable H2 and H3") +
  theme_usdt() +
  panel_gap
save_figure(fig_boundary, "Supplemental/Figure_S1.png", width = 6.9, height = 2.4)

# The same for the hypotheses, one row per parameter and one column per
# reliability (axes shared within each row)
reliab_show <- c(.1, .3, .5, .7)
hyp_labels  <- c(diff = "mu[I]-mu[A]~(H1)", rho = "rho~(H2)",
                 intercept = "beta[0]~(H3)", slope = "beta[1]~(H3)")
df_reliab <- df_true_est |>
  filter(parameter %in% names(hyp_labels), reliab %in% reliab_show) |>
  mutate(R_lab   = paste0("italic(R)=='", sub("^0\\.", ".", reliab), "'"),
         par_lab = hyp_labels[parameter])
ccc_reliab_labels <- ccc_reliab |>
  filter(parameter %in% names(hyp_labels), reliab %in% reliab_show) |>
  mutate(R_lab   = paste0("italic(R)=='", sub("^0\\.", ".", reliab), "'"),
         par_lab = hyp_labels[parameter],
         label   = paste("CCC =", sub("^0\\.", ".", sprintf("%.2f", CCC))))
fig_rows <- lapply(setNames(nm = names(hyp_labels)), function(p) {
  fig_p <- df_reliab |>
    filter(parameter == p) |>
    ggplot(aes(x = truth, y = mean_est)) +
    geom_abline(slope = 1, intercept = 0, linetype = "22", colour = "grey30",
                linewidth = 0.4) +
    geom_blank(aes(x = mean_est, y = truth)) +
    geom_point(size = 0.6, colour = "#256abf", alpha = 0.5) +
    geom_text(data = filter(ccc_reliab_labels, parameter == p),
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.08, vjust = 1.4, size = 2.3, inherit.aes = FALSE) +
    facet_grid(par_lab ~ R_lab, labeller = label_parsed) +
    labs(x = "True value", y = "Mean estimate") +
    theme_usdt() +
    theme(strip.text.y = element_text(angle = 0, hjust = 0))
  # Reliability labels only on the first row
  if (p != "diff") fig_p <- fig_p + theme(strip.text.x = element_blank())
  fig_p
})
fig_recovery_reliab <- patchwork::wrap_plots(fig_rows, ncol = 1) +
  patchwork::plot_layout(axis_titles = "collect")
save_figure(fig_recovery_reliab, "Supplemental/Figure_S2.png", width = 6.9, height = 6.2)

# Labels of the parameters
par_labels <- c(mu_dI     = "Influence~mean~mu[I]",
                diff      = "Difference~mu[I]-mu[A]~(H1)",
                sd_dD     = "Awareness~SD~sigma[A]",
                sd_dI     = "Influence~SD~sigma[I]",
                cov_d     = "Covariance~sigma[AI]",
                rho       = "Correlation~rho~(H2)",
                intercept = "Intercept~beta[0]~(H3)",
                slope     = "Slope~beta[1]~(H3)")

# Reference lines and Bradley's interval of the bias and coverage panels
measure_labels <- c(rel_bias = "Relative~bias", bias = "Bias",
                    coverage = "Coverage~of~the~'95%'~CI")
ref_lines <- data.frame(measure_lab = factor(measure_labels, levels = measure_labels),
                        y = c(0, 0, .95))
ref_band  <- data.frame(measure_lab = factor(measure_labels["coverage"], levels = measure_labels),
                        ymin = .925, ymax = .975)

# Relative bias and coverage of the random-effect SDs
fig_sd <- df_long |>
  filter(parameter %in% c("sd_dD", "sd_dI"), measure %in% c("rel_bias", "coverage")) |>
  summarise(value = mean(value), .by = c(parameter, measure, N, reliab)) |>
  mutate(par_lab     = par_labels[parameter],
         measure_lab = factor(measure_labels[measure], levels = measure_labels)) |>
  ggplot(aes(x = reliab, y = value, colour = factor(N))) +
  geom_rect(data = ref_band, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "grey92", inherit.aes = FALSE) +
  geom_hline(data = ref_lines[-2, ], aes(yintercept = y), linetype = "22",
             colour = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_grid(measure_lab ~ par_lab, scales = "free_y", switch = "y",
             labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  reliab_axis +
  labs(x = "Reliability", y = NULL) +
  theme_usdt() +
  theme(strip.placement = "outside")
save_figure(fig_sd, "Supplemental/Figure_S3.png", width = 5.2, height = 3.8)

# Bias and coverage of the correlation, per true correlation
fig_rho <- df_long |>
  filter(parameter == "rho", measure %in% c("bias", "coverage")) |>
  summarise(value = mean(value), .by = c(measure, N, reliab, rho)) |>
  mutate(rho_lab     = reorder(paste0("rho=='", sub("^0\\.", ".", rho), "'"), rho),
         measure_lab = factor(measure_labels[measure], levels = measure_labels)) |>
  ggplot(aes(x = reliab, y = value, colour = factor(N))) +
  geom_rect(data = ref_band, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            fill = "grey92", inherit.aes = FALSE) +
  geom_hline(data = ref_lines[-1, ], aes(yintercept = y), linetype = "22",
             colour = "grey30", linewidth = 0.5) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  facet_grid(measure_lab ~ rho_lab, scales = "free_y", switch = "y",
             labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  reliab_axis +
  labs(x = "Reliability", y = NULL) +
  theme_usdt() +
  theme(strip.placement = "outside") +
  panel_gap
save_figure(fig_rho, "Supplemental/Figure_S4.png", width = 6.9, height = 3.6)

# Power of H1 against the true difference, per trials per condition
fig_H1 <- df_long |>
  filter(parameter == "diff", measure %in% c("truth", "rejection")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  mutate(K_lab = paste(trials, "trials per condition")) |>
  ggplot(aes(x = truth, y = rejection, colour = factor(N))) +
  geom_hline(yintercept = .05, linetype = "22", colour = "grey30",
             linewidth = 0.5) +
  geom_point(size = 0.7) +
  facet_wrap(~ K_lab, nrow = 1) +
  scale_colour_manual(values = N_colours, name = "N") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = expression("True difference " * mu[I] - mu[A]),
       y = "Rejection rate of H1") +
  theme_usdt() +
  panel_gap
save_figure(fig_H1, "Supplemental/Figure_S5.png", width = 6.9, height = 2.6)

# Power of H1 (always testable) by trials and sample size for each true
# group-level difference (one per intercept and correlation), averaged over
# reliability, with the panels ordered by the absolute difference
# (supplemental material). Trials, not reliability, on the x axis: with fixed
# trials, higher reliability means larger true SDs, which make the
# group-level mean more variable
df_power_H1 <- df_long |>
  filter(parameter == "diff", measure %in% c("truth", "rejection_all")) |>
  pivot_wider(names_from = measure, values_from = value) |>
  summarise(truth = mean(truth), power = mean(rejection_all),
            .by = c(intercept, rho, N, trials)) |>
  mutate(truth = mean(truth), .by = c(intercept, rho)) |>
  arrange(abs(truth)) |>
  mutate(panel = paste0("mu[I]-mu[A]=='",
                        sub("^-", "−", sub("^(-?)0\\.", "\\1.", sprintf("%.3f", truth))), "'"),
         panel = factor(panel, levels = unique(panel)),
         label = ifelse(power > .995, "1", sub("^0\\.", ".", sprintf("%.2f", power))))
fig_power_H1 <- ggplot(df_power_H1, aes(x = factor(trials), y = factor(N), fill = power)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = label, colour = power > .55,
                fontface = ifelse(power >= .80, "bold", "plain")), size = 2.1) +
  facet_wrap(~ panel, ncol = 5, labeller = label_parsed) +
  scale_fill_gradientn(colours = c("#eef4fc", "#86b6ef", "#256abf", "#0d366b"),
                       limits = c(0, 1), breaks = c(0, .5, .8, 1),
                       labels = c("0", ".50", ".80", "1"), name = "Power") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  labs(x = "Trials per condition (K)", y = "Sample size (N)") +
  theme_usdt(base_size = 8) +
  theme(panel.border = element_blank(), axis.line = element_blank(),
        axis.ticks = element_blank(),
        legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.25, "cm"),
        plot.margin = margin(5.5, 5.5, 14, 5.5))
save_figure(fig_power_H1, "Supplemental/Figure_S6.png", width = 6.9, height = 4.6)

# Rejection rates of H2 and H3 over the testable fits and over all replications
rejection_lines <- scale_linetype_manual(
  values = c(rejection = "solid", rejection_all = "22"),
  labels = c("Testable fits", "All replications"), name = NULL)
fig_H2 <- df_long |>
  filter(parameter == "rho", measure %in% c("rejection", "rejection_all")) |>
  summarise(value = mean(value), .by = c(measure, N, reliab, rho)) |>
  mutate(rho_lab = reorder(paste0("rho=='", sub("^0\\.", ".", rho), "'"), rho)) |>
  ggplot(aes(x = reliab, y = value, colour = factor(N), linetype = measure)) +
  geom_hline(yintercept = .05, colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ rho_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  rejection_lines +
  reliab_axis +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = NULL, y = "Rejection rate of H2") +
  theme_usdt() +
  panel_gap
fig_H3 <- df_long |>
  filter(parameter == "intercept", measure %in% c("rejection", "rejection_all")) |>
  summarise(value = mean(value), .by = c(measure, N, reliab, intercept)) |>
  mutate(a_lab = reorder(paste0("beta[0]=='", sub("^0\\.", ".", intercept), "'"), intercept)) |>
  ggplot(aes(x = reliab, y = value, colour = factor(N), linetype = measure)) +
  geom_hline(yintercept = .05, colour = "grey60", linewidth = 0.4) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ a_lab, nrow = 1, labeller = label_parsed) +
  scale_colour_manual(values = N_colours, name = "N") +
  rejection_lines +
  reliab_axis +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(x = "Reliability", y = "Rejection rate of H3") +
  theme_usdt() +
  panel_gap
fig_H2H3 <- patchwork::wrap_plots(A = fig_H2, B = fig_H3,
                                  design = "AAAAA\nBBB##", guides = "collect")
save_figure(fig_H2H3, "Supplemental/Figure_S7.png", width = 6.9, height = 4.4)

# Power over the testable fits, made with Figure 7 in Section 10
save_figure(fig_power$testable, "Supplemental/Figure_S8.png", width = 6.9, height = 4.4)

# ─────────────────────────────────────────────────────────────────────────────
