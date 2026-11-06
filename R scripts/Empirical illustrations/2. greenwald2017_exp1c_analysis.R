# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : greenwald2017_exp1c.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 02-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Greenwald and De Houwer (2017), Experiment 1c (data set S18).
# Data: https://osf.io/ydeza/.
#   1. Awareness task: CS visibility test (144 trials).
#   2. Influence task: conditioning test (96 trials).
#
# The trial-level data and the authors' d' come from
# "R scripts/Empirical ilustrations/1. greenwald2017_prepare_data.R", which must
# be run first. Paths are relative to the project root.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load packages and user-defined functions
# ─────────────────────────────────────────────────────────────────────────────

library(uSDT)
library(dplyr)

# Errors-in-variables R functions (Klauer et al., 1998)
# Retrieved from: https://osf.io/36hk2/files/osfstorage (in the .zip folder!)
source("R scripts/R functions/eiv_functions.R")

# Shared helper to save figures in Figures/
source("R scripts/R functions/plot_functions.R")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Load data
# ─────────────────────────────────────────────────────────────────────────────

# Load data and observed dprimes from greenwald
greenwald2017 <- read.csv("Data/Greenwald 2017/greenwald2017_long.csv")
dprimes       <- read.csv("Data/Greenwald 2017/greenwald2017_dprimes.csv")

# The authors' d' of the participants in their regression (N = 350)
s18 <- dprimes |>
  filter(study == "S18", !is.na(direct), !is.na(indirect))

# Analysed trials of the participants with both tasks. The signal is the CS
# paired with the right key, and a right-key response is a "yes".
trials <- greenwald2017 |>
  filter(study == "S18", exclusion == "none", !is.na(response)) |>
  filter(n_distinct(task) == 2, .by = subject)

stopifnot(setequal(s18$subject, unique(trials$subject))) # Small check

# Prepare uSDT data
d <- usdt_data_long(trials,
                    task_col         = "task",
                    task_levels      = c(awareness = "direct", influence = "indirect"),
                    subject_col      = "subject",
                    condition_col    = "cs",
                    condition_levels = c(signal = "R", noise = "L"),
                    response_col     = "response",
                    response_levels  = c(signal = "R", noise = "L"))

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Standard two-step analysis without accounting for meas. error
# ─────────────────────────────────────────────────────────────────────────────

# Mean d' of each task
colMeans(s18[c("direct", "indirect")])

# Hypothesis 1: group-level difference (Awareness - Influence)
t.test(x = s18[,"direct"], y = s18[,"indirect"], paired = TRUE)

# Hypothesis 2: correlation between d-primes
cor.test(x = s18[,"direct"], y = s18[,"indirect"])

# Hypothesis 3: linear regression between d-primes
ols <- lm(indirect ~ direct, data = s18)
summary(ols)
confint(ols, level = .99)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Errors-in-variables regression
# ─────────────────────────────────────────────────────────────────────────────

# Set seed for starting values in eiv
set.seed(2017)

# Estimate regression intercept
eiv_fit <- eiv_lm(s18$direct, s18$indirect)

# Results on intercept and slope
eiv <- eiv_fit |>
  filter(parameter %in% c("intercept", "slope")) |>
  mutate(conf.low  = estimate - qt(.975, nrow(s18) - 2) * se,
         conf.high = estimate + qt(.975, nrow(s18) - 2) * se,
         p.value   = round(
           2 * pt(-abs(estimate / se), df = nrow(s18) - 2),3))
print(eiv, digits = 3)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Multilevel SDT model, frequentist estimation
# ─────────────────────────────────────────────────────────────────────────────

# The fit reaches the boundary: the SD of the indirect d' is estimated at
# zero, so rho and the intercept have no standard errors
m <- hsdt(d)
summary(m)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Multilevel SDT model, Bayesian estimation
# ─────────────────────────────────────────────────────────────────────────────

# We define priors using latent regression parameterization
priors <- usdt_priors(intercept = "normal(0, 1)", cor_dprime = "lkj(1)")
priors

# Sampling from the model (more draws for stable bayes factors)
m_bayes <- hsdt(
  data = d,
  estimation = "bayesian",
  latent = "regression",
  priors = priors,
  backend = "cmdstanr",
  chains = 4,
  iter = 27000,
  warmup = 2000,
  seed = 2017)

# Save fitted model
saveRDS(m_bayes, file = "Fits/bhsdt_greenwald2017.rds")

# Model summary with hypotesis H1 to H3.
summary(m_bayes)

# Bayes factors via Savage-Dickey density ratio test
# Hypothesis 1: group-level difference (Influence - Awareness)
usdt_bf(fit = m_bayes, hypothesis = "diff = 0")

# Hypothesis 2: correlation between d-primes
usdt_bf(fit = m_bayes, hypothesis = "rho = 0")

# Hypothesis 3: linear regression between d-primes
usdt_bf(fit = m_bayes, hypothesis = "intercept = 0")

# Save all results in one objetct
set.seed(2017)
bf <- usdt_bf(m_bayes, plot = FALSE)

# Reliability of the observed d' in each task (and participant!)
reliability <- usdt_reliability(m_bayes)
reliability

# Posterior mean and 95% CrI of the mean d' of each task
mu_d <- as.data.frame(m_bayes$fit,
                      variable = c("b_dA_Intercept", "b_dI_Intercept")) |>
  sapply(function(x) c(estimate = mean(x), quantile(x, c(.025, .975))))
colnames(mu_d) <- c("direct", "indirect")

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Figures for the main text
# ─────────────────────────────────────────────────────────────────────────────

# Population ROC curves of each task (frequentist fit), with 95% Wald bands.
# The captions are dropped because the manuscript gives its own notes.
fig_roc <- plot(m, type = "roc", caption = FALSE)
save_figure(fig_roc, "Manuscript/Figure_8.png", width = 5, height = 5)

# Observed vs. model-estimated sensitivities (Bayesian fit). Panels are stacked
# and use the full page width, so that the shrinkage is easy to see.
fig_reg <- plot(m_bayes, type = "regression", nrow = 2, caption = FALSE)
save_figure(fig_reg, "Manuscript/Figure_9.png", width = 6.9, height = 7.5)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Save results
# ─────────────────────────────────────────────────────────────────────────────

# This part of the code save the results used in the main text. We don't save
# model fits due to their size.
dir.create("Results/Empirical illustrations", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n           = nrow(s18),
             ols         = data.frame(
               term      = c("intercept", "slope"),
               estimate  = unname(coef(ols)),
               conf.low  = confint(ols, level = .99)[, 1],
               conf.high = confint(ols, level = .99)[, 2],
               level     = .99),
             eiv         = eiv,
             ml          = m$tests,
             bayes       = m_bayes$tests,
             bf          = data.frame(bf)[c("hypothesis", "estimate", "conf.low",
                                            "conf.high", "BF10", "BF01")],
             mu_d        = mu_d,
             reliability = reliability$tasks,
             diagnostics = m_bayes$diagnostics),
        "Results/Empirical illustrations/greenwald2017_exp1c.rds")

# ─────────────────────────────────────────────────────────────────────────────
