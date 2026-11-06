# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : vadillo2025_exp2.R                                         ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 18-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Vadillo, Malejka, and Shanks (2025), Experiment 2.
# Data: https://osf.io/jp3gx/. Also included in {uSDT} package
#   1. Influence task: visual search (384 trials).
#   2. Awareness task: Old/New recognition test (64 trials).
#
# Trial exclusions and d' follow the authors' scripts (02_initial_analyses.R,
# filtering_functions.R, and scoring_functions.R on the same OSF project).

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load packages
# ─────────────────────────────────────────────────────────────────────────────

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
data_dir <- "Data/Vadillo 2024"
osf_url  <- c(cuing.csv = "https://osf.io/download/6b2fe/",
              aware.csv = "https://osf.io/download/bpxe6/")
osf_md5  <- c(cuing.csv = "f549cf261cf80570f1bc3df1a6f8d952",
              aware.csv = "f62b49d5b8d043a4e998538e307c289c")

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

# Load both dataframes
cuing <- read.csv(file.path(data_dir, "cuing.csv"))
aware <- read.csv(file.path(data_dir, "aware.csv"))

# Valid search trials: correct responses, RTs up to 10 s, and RTs within 3 SDs
# of the participant's mean of the remaining trials
cuing <- cuing |>
  mutate(valid = acc == 1 & rt <= 10000) |>
  mutate(m_rt  = mean(rt[valid]),
         sd_rt = sd(rt[valid]), .by = subj) |>
  mutate(valid = valid & rt >= m_rt - 3 * sd_rt & rt <= m_rt + 3 * sd_rt) |>
  select(-m_rt, -sd_rt)

# uSDT data. Old displays are the signal condition, and RTs are split at each
# participant's median.
d <- usdt_data_tasks(
  awareness           = aware,
  influence         = filter(cuing, valid),
  subject_col      = "subj",
  condition_col    = "condition",
  condition_levels = c(signal = "old", noise = "new"),
  response_col     = list(awareness = "judged.old", influence = "rt"),
  response_levels  = list(awareness   = c(signal = 1, noise = 0),
                          influence = c(signal = "faster", noise = "slower")),
  dichotomize      = list(awareness = FALSE, influence = TRUE)
)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Original analysis
# ─────────────────────────────────────────────────────────────────────────────

# Contextual cuing in ms (new minus old mean RT) against zero
cuing_ms <- cuing |>
  filter(valid) |>
  summarise(rt = mean(rt), .by = c(subj, condition)) |>
  pivot_wider(names_from = condition, values_from = rt) |>
  mutate(effect = new - old)

t_cuing <- t.test(cuing_ms$effect)
t_cuing

# Recognition d' against zero, with the log-linear correction (32 trials per
# condition)
recognition <- aware |>
  summarise(judged.old = mean(judged.old), .by = c(subj, condition)) |>
  pivot_wider(names_from = condition, values_from = judged.old) |>
  mutate(dprime = qnorm((old * 32 + .5) / 33) - qnorm((new * 32 + .5) / 33))

t_recognition <- t.test(recognition$dprime)
t_recognition

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: The three hypotheses with observed d' (two-step approach)
# ─────────────────────────────────────────────────────────────────────────────

# Compute observed d-primes
scores <- sdt_moments(d) |>
  select(subj, task, dprime) |>
  pivot_wider(names_from = task, values_from = dprime)

# Hypothesis 1: group-level differences (Influence - Awareness)
h1 <- t.test(scores$Influence, scores$Awareness, paired = TRUE)
h1

# Hypothesis 2: correlation between d-primes
h2 <- cor.test(scores$Awareness, scores$Influence)
h2

# Hypothesis 3: linear regression between d-primes
h3 <- lm(Influence ~ Awareness, data = scores)
summary(h3)

# Save two-step approach results
two_step <- data.frame(
  hypothesis = c("H1", "H2", "H3"),
  estimate   = c(unname(h1$estimate), unname(h2$estimate), coef(h3)[[1]]),
  conf.low   = c(h1$conf.int[1], h2$conf.int[1], confint(h3)[1, 1]),
  conf.high  = c(h1$conf.int[2], h2$conf.int[2], confint(h3)[1, 2]),
  statistic  = c(h1$statistic, h2$statistic, summary(h3)$coefficients[1, 3]),
  df         = c(h1$parameter, h2$parameter, h3$df.residual),
  p.value    = c(h1$p.value, h2$p.value, summary(h3)$coefficients[1, 4])
)
two_step

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
  seed = 2025,
  control = list(adapt_delta = .99))

# Save fitted model
saveRDS(m_bayes, file = "Fits/bhsdt_vadillo2025.rds")

# Model summary with hypotesis H1 to H3.
summary(m_bayes)

# Bayes factors via Savage-Dickey density ratio test
# Hypothesis 1: group-level difference (Awareness - Influence)
usdt_bf(fit = m_bayes, hypothesis = "diff = 0")

# Hypothesis 2: correlation between d-primes
usdt_bf(fit = m_bayes, hypothesis = "rho = 0")

# Hypothesis 3: linear regression between d-primes
usdt_bf(fit = m_bayes, hypothesis = "intercept = 0")

# Save all results in one objetct
set.seed(2025)
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

# Observed vs. model-estimated sensitivities (Bayesian fit). Panels are stacked
# and use the full page width; the caption is dropped because the manuscript
# gives its own note.
fig_reg <- plot(m_bayes, type = "regression", nrow = 2, caption = FALSE)
save_figure(fig_reg, "Manuscript/Figure_10.png", width = 6.9, height = 7.5)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 9: Save results
# ─────────────────────────────────────────────────────────────────────────────

# This part of the code save the results used in the main text. We don't save
# model fits due to their size.
dir.create("Results/Empirical illustrations", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(n           = nrow(scores),
             excluded    = mean(!cuing$valid),
             cuing_ms    = c(estimate = unname(t_cuing$estimate),
                             t = unname(t_cuing$statistic),
                             df = unname(t_cuing$parameter),
                             p.value = t_cuing$p.value),
             recognition = c(estimate = unname(t_recognition$estimate),
                             t = unname(t_recognition$statistic),
                             df = unname(t_recognition$parameter),
                             p.value = t_recognition$p.value),
             two_step    = two_step,
             ml          = m$tests,
             bayes       = m_bayes$tests,
             bf          = data.frame(bf)[c("hypothesis", "estimate", "conf.low",
                                            "conf.high", "BF10", "BF01")],
             mu_d        = mu_d,
             reliability = reliability$tasks,
             diagnostics = m_bayes$diagnostics),
        "Results/Empirical illustrations/vadillo2025_exp2.rds")

# ─────────────────────────────────────────────────────────────────────────────
