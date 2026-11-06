# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : wentura2025_exp2a.R                                        ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 18-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Wentura, Rohr, and Kiefer (2025), Experiment 2a.
# Data: https://osf.io/9t4pk/ (folder "Experiments 2ab"). 360 trials per task
#   1. Influence task: emotion categorization of a visible target face.
#   2. Awareness task: emotion categorization of the masked prime.
#
# Paths are relative to the project root.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load packages
# ─────────────────────────────────────────────────────────────────────────────

# pak::pak("RicardoReySaez/uSDT@dev")
library(uSDT)
library(dplyr)
library(tidyr)

# Shared helper to save figures in Figures/
source("R scripts/R functions/plot_functions.R")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Download and prepare data
# ─────────────────────────────────────────────────────────────────────────────

# Data directory and osf urls
data_dir <- "Data/Wentura 2025"
osf_url  <- c(longformat_Exp2a.csv           = "https://osf.io/download/u47c2/",
              longformat_awareness_Exp2a.csv = "https://osf.io/download/gpzbd/")
osf_md5  <- c(longformat_Exp2a.csv           = "40dfa7e86a81ef4add6526350035ae62",
              longformat_awareness_Exp2a.csv = "39dcc9b6abe89f5ea837ae884dcfe613")

# Create directory and save osf data
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
for (file in names(osf_url)) {
  path <- file.path(data_dir, file)
  if (!file.exists(path) || unname(tools::md5sum(path)) != osf_md5[[file]]) {
    download.file(osf_url[[file]], path, mode = "wb")
  }
  if (unname(tools::md5sum(path)) != osf_md5[[file]]) {
    stop("The file downloaded from OSF is not the ", file, " analysed here ",
         "(MD5 ", osf_md5[[file]], ").", call. = FALSE)
  }
}

# Files are UTF-8 with a byte-order mark, and missing values are blank spaces
read_wentura <- function(file) {
  read.csv(file.path(data_dir, file),
           fileEncoding = "UTF-8-BOM", na.strings = c("", " "))
}

priming   <- read_wentura("longformat_Exp2a.csv")
awareness <- read_wentura("longformat_awareness_Exp2a.csv")

# Indirect task. rt4 is the authors' RT, with errors and outliers set to
# missing. Congruent trials (c = 1) are the signal condition.
rt_df <- priming |>
  filter(!is.na(rt4)) |>
  select(subj = vp, cond = c, rt = rt4)

# Share of trials discarded by the authors
excluded <- mean(is.na(priming$rt4))
excluded

# Direct task. Sad primes are the signal condition, and the response is 1 when
# the participant answered "sad".
aw_df <- awareness |>
  mutate(cond = as.integer(pe == "sad"),
         resp = ifelse(pe == "sad", hit, fa)) |>
  select(subj = vp, cond, resp)

# uSDT data. RTs are split at each participant's median.
d <- usdt_data_tasks(
  awareness        = aw_df,
  influence        = rt_df,
  subject_col      = "subj",
  condition_col    = "cond",
  condition_levels = c(signal = 1, noise = 0),
  response_col     = list(awareness = "resp", influence = "rt"),
  response_levels  = list(awareness   = c(signal = 1, noise = 0),
                          influence = c(signal = "faster", noise = "slower")),
  dichotomize      = list(awareness = FALSE, influence = TRUE)
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Original analysis
# ─────────────────────────────────────────────────────────────────────────────

# Priming effect in ms (incongruent minus congruent mean RT)
priming_ms <- rt_df |>
  summarise(rt = mean(rt), .by = c(subj, cond)) |>
  pivot_wider(names_from = cond, values_from = rt, names_prefix = "rt_") |>
  mutate(effect = rt_0 - rt_1)

# Observed d' of each task, joined to the priming effect by participant
scores <- sdt_moments(d) |>
  select(subj, task, dprime) |>
  pivot_wider(names_from = task, values_from = dprime) |>
  mutate(subj = as.character(subj)) |>
  left_join(priming_ms |> transmute(subj = as.character(subj), effect),
            by = "subj")

# Effect sizes (dz)
mean(scores$effect) / sd(scores$effect)
mean(scores$Awareness) / sd(scores$Awareness)

# Priming effect and direct d' against zero, and their correlation
t_priming <- t.test(scores$effect)
t_awareness  <- t.test(scores$Awareness)
r_priming <- cor.test(scores$effect, scores$Awareness)
t_priming; t_awareness; r_priming

# 95% CIs of the effect sizes reported by the authors (N = 86)
ci_dz <- function(dz, n, level = .95) {
  t_obs <- dz * sqrt(n)
  df <- n - 1
  a <- (1 - level) / 2

  cdf <- function(ncp) {
    integrate(function(v) {
      pnorm(t_obs * sqrt(v / df) - ncp) * dchisq(v, df)
    }, 0, Inf, rel.tol = 1e-9)$value
  }

  ncp <- sapply(c(1 - a, a), function(p) {
    uniroot(function(ncp) cdf(ncp) - p,
            c(t_obs - 10, t_obs + 10), tol = 1e-8)$root
  })

  ncp / sqrt(n)
}

# Also 95% CIs for correlations
ci_r <- function(r, n, level = .95) {
  tanh(atanh(r) + c(-1, 1) * qnorm(1 - (1 - level) / 2) / sqrt(n - 3))
}

# Reported CIs
reported_ci <- rbind(priming     = ci_dz(0.34, 86),
                     direct      = ci_dz(0.27, 86),
                     correlation = ci_r(-0.03, 86))
round(reported_ci, 2)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: The three hypotheses with observed d' (two-step approach)
# ─────────────────────────────────────────────────────────────────────────────

# Hypothesis 1: group-level differences
t.test(scores$Influence, scores$Awareness, paired = TRUE)

# Hypothesis 2: correlation
r_dprime <- cor.test(scores$Awareness, scores$Influence)
r_dprime

# Hypothesis 3: intercept
ols <- lm(Influence ~ Awareness, data = scores)
summary(ols)
confint(ols)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Multilevel SDT model, frequentist estimation
# ─────────────────────────────────────────────────────────────────────────────

# Estimate hierarchical SDT model
m <- hsdt(d)

# Summary and hypothesis H1 to H3
summary(m)

# Reliability estimates
usdt_reliability(m)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Multilevel SDT model, Bayesian estimation
# ─────────────────────────────────────────────────────────────────────────────

# Priors that mirror the original Bayes factors: JSZ, a Cauchy(0, r * SD)
prior_scale <- function(r) {
  round(r * c(awareness = sd(scores$Awareness),
              influence = sd(scores$Influence)), 3)
}

# Function to automatically set the priors for each measure in uSDT format
matched_priors <- function(r) {
  s <- prior_scale(r)
  usdt_priors(dprime = list(
    awareness = sprintf("cauchy(0, %.3f)", s[["awareness"]]),
    influence = sprintf("cauchy(0, %.3f)", s[["influence"]])))
}

# The equivalent to JSZ prior
r0 <- sqrt(2) / 2
matched_priors(r0)

# Sampling from the model (more draws for stable bayes factors)
m_bayes <- hsdt(
  data = d,
  estimation = "bayesian",
  priors = matched_priors(r0),
  backend = "cmdstanr",
  warmup = 2000,
  iter = 27000,
  seed = 2025)

# Save fitted model
saveRDS(m_bayes, file = "Fits/bhsdt_wentura2025.rds")

# Model summary with hypotesis H1 to H3.
summary(m_bayes)

# Bayes factors via Savage-Dickey density ratio test
# Hypothesis 1: group-level difference (Awareness - Influence)
usdt_bf(fit = m_bayes, hypothesis = "diff = 0")

# Hypothesis 2: correlation between d-primes
usdt_bf(fit = m_bayes, hypothesis = c("rho = 0", "rho > 0"))

# Hypothesis 3: linear regression between d-primes
usdt_bf(fit = m_bayes, hypothesis = "intercept = 0")

# Save all results in one objetct
set.seed(2025)
bf <- usdt_bf(m_bayes, plot = FALSE)

# Reliability of the observed d' in each task (and participant!)
reliability <- usdt_reliability(m_bayes)
reliability

# Mean d' of each task: posterior mean, 95% CrI, and Bayes factor against
# zero (Savage-Dickey)
posterior_mu_d <- as.data.frame(
  m_bayes$fit,
  variable = c("b_dA_Intercept", "b_dI_Intercept"))
names(posterior_mu_d) <- c("awareness", "influence")

# Function to compute posterior density
post_dens_0  <- sapply(posterior_mu_d, function(x) {
  logspline::dlogspline(0, logspline::logspline(x))
})

# Prior density value for each measure at 0
prior_dens_0 <- dcauchy(0, 0, prior_scale(r0))

# Bayes factors on simple effects
BF01_mu      <- post_dens_0 / prior_dens_0

# Save results
mu_d <- rbind(sapply(posterior_mu_d, function(x)
  c(estimate = mean(x), quantile(x, c(.025, .975)))),
  BF01 = BF01_mu, BF10 = 1 / BF01_mu)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Sensitivity analysis to the prior scale
# ─────────────────────────────────────────────────────────────────────────────

# The three values of r used by the autors
r_vals <- c(sqrt(2) / 2, 1, sqrt(2))

# Re-estimate the model changing only these prior values
sensitivity <- lapply(r_vals, function(r) {
  # Refit the model with new priors
  fit <- hsdt(
    data = d,
    estimation = "bayesian",
    priors = matched_priors(r),
    backend = "cmdstanr",
    warmup = 2000,
    iter = 27000,
    seed = 2025)

  # Compute Bayes Factors
  set.seed(2025)
  usdt_bf(fit, plot = FALSE) |>
    transmute(r = round(r, 3), hypothesis, estimate, conf.low, conf.high,
              BF10, BF01)
})

# Save final sensitivity analysis results
sensitivity <- bind_rows(sensitivity)
sensitivity |>
  mutate(across(where(is.numeric), \(x) round(x, 2)))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Figures for the main text
# ─────────────────────────────────────────────────────────────────────────────

# Prior and posterior of rho for the point-null (Savage-Dickey) and the
# directional hypothesis, side by side. The caption is dropped because the
# manuscript gives its own note.
set.seed(2025)
bf_rho  <- usdt_bf(m_bayes, hypothesis = c("rho = 0", "rho > 0"), plot = FALSE)
fig_rho <- plot(bf_rho, caption = FALSE)
save_figure(fig_rho, "Manuscript/Figure_11.png", width = 7.2, height = 4.1)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 10: Save results
# ─────────────────────────────────────────────────────────────────────────────

# Function to store test values
test_values <- function(x) {
  c(estimate = unname(x$estimate), statistic = unname(x$statistic),
    df = unname(x$parameter), p.value = x$p.value)
}

# Save final results
dir.create("Results/Empirical illustrations", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n           = nrow(scores),
             excluded    = excluded,
             priming     = test_values(t_priming),
             awareness   = test_values(t_awareness),
             r_priming   = test_values(r_priming),
             r_dprime    = test_values(r_dprime),
             ols         = data.frame(term      = c("intercept", "slope"),
                                      estimate  = unname(coef(ols)),
                                      conf.low  = confint(ols)[, 1],
                                      conf.high = confint(ols)[, 2],
                                      p.value   = summary(ols)$coefficients[, 4],
                                      level     = .95),
             reported_ci = reported_ci,
             sd_dprime   = c(awareness = sd(scores$Awareness),
                             influence = sd(scores$Influence)),
             prior_scale = prior_scale(r0),
             bayes       = m_bayes$tests,
             bf          = data.frame(bf)[c("hypothesis", "estimate", "conf.low",
                                            "conf.high", "BF10", "BF01")],
             mu_d        = mu_d,
             reliability = reliability$tasks,
             sensitivity = data.frame(sensitivity),
             diagnostics = m_bayes$diagnostics),
        "Results/Empirical illustrations/wentura2025_exp2a.rds")

# ─────────────────────────────────────────────────────────────────────────────
