# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : greenwald2017_series1.R                                    ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 06-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Greenwald and De Houwer (2017), all 18 data sets of Series 1 (S1 to S18).
# Data: https://osf.io/ydeza/.
#   1. Awareness task: CS visibility test.
#   2. Influence task: conditioning test.
#
# The three hypotheses in each data set, with the same analyses that the main
# text reports for S18 (Experiment 1c) in
# "R scripts/Empirical illustrations/2. greenwald2017_exp1c_analysis.R":
#   - the two-step approach with the authors' d' (Section 4),
#   - the errors-in-variables regression of Klauer et al. (1998) on the
#     authors' d' (Section 5),
#   - the multilevel SDT model by maximum likelihood (Section 6), and
#   - the multilevel SDT model by Bayesian estimation, with the latent
#     regression parameterization and the same priors (Section 7).

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load packages and user-defined functions
# ─────────────────────────────────────────────────────────────────────────────

library(uSDT)
library(dplyr)
library(tidyr)

# Klauer et al. (1998) EIV functions: https://osf.io/36hk2/
source("R scripts/R functions/eiv_functions.R")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Load data
# ─────────────────────────────────────────────────────────────────────────────

# Trial data and authors' observed d' for all 18 studies.
greenwald2017 <- read.csv("Data/Greenwald 2017/greenwald2017_long.csv")
dprimes       <- read.csv("Data/Greenwald 2017/greenwald2017_dprimes.csv") |>
  filter(!is.na(direct), !is.na(indirect))
studies       <- paste0("S", 1:18)

# Match participants with both tasks to the authors' d' sample.
usdt_study <- function(study_id) {
  trials <- greenwald2017 |>
    filter(study == study_id, exclusion == "none", !is.na(response)) |>
    filter(n_distinct(task) == 2, .by = subject)
  stopifnot(setequal(unique(trials$subject),
                     dprimes$subject[dprimes$study == study_id]))
  usdt_data_long(trials,
                 task_col         = "task",
                 task_levels      = c(awareness = "direct", influence = "indirect"),
                 subject_col      = "subject",
                 condition_col    = "cs",
                 condition_levels = c(signal = "R", noise = "L"),
                 response_col     = "response",
                 response_levels  = c(signal = "R", noise = "L"))
}

d <- setNames(lapply(studies, usdt_study), studies)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: The three hypotheses with observed d' (two-step approach)
# ─────────────────────────────────────────────────────────────────────────────

# Observed d': H1 paired difference; H2 Pearson r; H3 intercept.
two_step <- function(study_id) {
  obs <- filter(dprimes, study == study_id)
  h1  <- t.test(obs$indirect, obs$direct, paired = TRUE)
  h2  <- cor.test(obs$direct, obs$indirect)
  h3  <- lm(indirect ~ direct, data = obs)
  data.frame(
    n          = nrow(obs),
    hypothesis = c("H1", "H2", "H3"),
    estimate   = c(unname(h1$estimate), unname(h2$estimate), coef(h3)[[1]]),
    conf.low   = c(h1$conf.int[1], h2$conf.int[1], confint(h3)[1, 1]),
    conf.high  = c(h1$conf.int[2], h2$conf.int[2], confint(h3)[1, 2]),
    p.value    = c(h1$p.value, h2$p.value, summary(h3)$coefficients[1, 4])
  )
}

observed <- lapply(setNames(studies, studies), two_step) |>
  bind_rows(.id = "study") |>
  mutate(approach = "two-step")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Errors-in-variables regression
# ─────────────────────────────────────────────────────────────────────────────

# H3 only; t-based 95% CIs and p values (N - 2 df).
# Seed reproduces S18; S1/S4/S14 may have boundary fits (no SE).
eiv <- lapply(setNames(studies, studies), function(study_id) {
  obs <- filter(dprimes, study == study_id)
  df  <- nrow(obs) - 2
  set.seed(2017)
  eiv_lm(obs$direct, obs$indirect) |>
    filter(parameter == "intercept") |>
    transmute(hypothesis = "H3", estimate,
              conf.low  = estimate - qt(.975, df) * se,
              conf.high = estimate + qt(.975, df) * se,
              p.value   = 2 * pt(-abs(estimate / se), df),
              note      = ifelse(is.na(se), "boundary", ""))
}) |>
  bind_rows(.id = "study") |>
  left_join(distinct(observed, study, n), by = "study") |>
  mutate(approach = "eiv")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Multilevel SDT model, frequentist estimation
# ─────────────────────────────────────────────────────────────────────────────

# Boundary fits may lack SEs for H2/H3; mark them in the tables.
fits <- lapply(d, function(data) suppressWarnings(hsdt(data)))

frequentist <- lapply(fits, function(m) {
  m$tests |>
    filter(term != "slope") |>
    transmute(hypothesis, estimate, conf.low, conf.high, p.value,
              note = case_when(status == "ok" ~ "",
                               grepl("boundary", reason) ~ "boundary",
                               .default = reason))
}) |>
  bind_rows(.id = "study") |>
  left_join(distinct(observed, study, n), by = "study") |>
  mutate(approach = "frequentist")

# Mean d' by task (A: awareness; I: influence).
round(t(sapply(fits, function(m) lme4::fixef(m$fit)[c("d_A", "d_I")])), 3)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Multilevel SDT model, Bayesian estimation
# ─────────────────────────────────────────────────────────────────────────────

# Latent regression; main-text priors for intercept and correlation.
# Remaining priors: uSDT defaults; posterior p values are two-sided.
priors <- usdt_priors(intercept = "normal(0, 1)", cor_dprime = "lkj(1)")
priors

# Four chains: 2,000 warmup + 2,500 draws each (vs. 25,000).
fits_bayes <- lapply(d, function(data) {
  hsdt(data, estimation = "bayesian", latent = "regression", priors = priors,
       backend = "cmdstanr", chains = 4, iter = 4500, warmup = 2000,
       seed = 2017, refresh = 0)
})

# Convergence: max R-hat, min bulk/tail ESS, and divergences.
lapply(fits_bayes, function(m) {
  with(m$diagnostics, data.frame(max_rhat     = round(max_rhat, 3),
                                 min_ess_bulk = round(min_ess_bulk),
                                 min_ess_tail = round(min_ess_tail),
                                 divergent    = divergent))
}) |>
  bind_rows(.id = "study")

# Savage–Dickey Bayes factors; seed fixes density estimates.
set.seed(2017)
bayesian <- lapply(fits_bayes, function(m) {
  m$tests |>
    filter(term != "slope") |>
    transmute(hypothesis, estimate, conf.low, conf.high, p.value,
              BF01 = usdt_bf(m, plot = FALSE)$BF01)
}) |>
  bind_rows(.id = "study") |>
  left_join(distinct(observed, study, n), by = "study") |>
  mutate(approach = "bayesian")

# Posterior means of task-level d' means and between-subject SDs.
pop_pars <- c(mu_A    = "b_dA_Intercept",        mu_I    = "b_dI_Intercept",
              sigma_A = "sd_subj__dA_Intercept", sigma_I = "sd_subj__dI_Intercept")
round(t(sapply(fits_bayes, function(m) {
  setNames(colMeans(as.data.frame(m$fit, variable = unname(pop_pars))),
           names(pop_pars))
})), 3)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Reliability of each measure
# ─────────────────────────────────────────────────────────────────────────────

# Reliability: true between-subject d' variance / total variance.
# ML: point estimates; Bayesian: posterior means and 95% intervals.
reliability <- bind_rows(
  frequentist = lapply(fits, function(m) usdt_reliability(m)$tasks) |>
    bind_rows(.id = "study"),
  bayesian    = lapply(fits_bayes, function(m) usdt_reliability(m)$tasks) |>
    bind_rows(.id = "study"),
  .id = "approach"
)

reliability |>
  mutate(study = factor(study, levels = studies),
         value = case_when(
           approach == "bayesian" ~ sprintf("%.3f [%.3f, %.3f]", reliability,
                                            conf.low, conf.high),
           .default = sprintf("%.3f", reliability))) |>
  select(study, task, approach, value) |>
  pivot_wider(names_from = c(task, approach), values_from = value) |>
  arrange(study) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Results
# ─────────────────────────────────────────────────────────────────────────────

# H1: d' difference; H2: correlation; H3: intercept.
# Estimates, 95% intervals, p values, and Bayesian BF01.
approaches <- c("two-step", "eiv", "frequentist", "bayesian")

results <- bind_rows(observed, eiv, frequentist, bayesian) |>
  mutate(study = factor(study, levels = studies),
         value = case_when(
           !is.na(conf.low) ~ sprintf("%.3f [%.3f, %.3f]", estimate,
                                      conf.low, conf.high),
           .default = sprintf("%.3f (%s)", estimate, note)),
         p     = ifelse(is.na(p.value), "",
                        format.pval(p.value, digits = 2, eps = .001)))

for (h in c("H1", "H2", "H3")) {
  cat("\n", h, "\n")
  results |>
    filter(hypothesis == h) |>
    select(study, n, approach, value, p) |>
    pivot_wider(names_from = approach, values_from = c(value, p),
                names_glue = "{approach}_{.value}") |>
    select(study, n, any_of(paste0(rep(approaches, each = 2), c("_value", "_p")))) |>
    arrange(study) |>
    left_join(filter(bayesian, hypothesis == h) |>
                transmute(study, BF01 = signif(BF01, 3)), by = "study") |>
    as.data.frame() |>
    print(row.names = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────
