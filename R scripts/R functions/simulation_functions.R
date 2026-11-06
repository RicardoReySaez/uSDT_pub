# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : simulation_functions.R                                     ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 22-09-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Helper functions of the validation and robustness simulation studies.
#
# The functions are grouped as follows:
#   Section 2: population parameters of the validation simulation study,
#   Section 3: population parameters and reliability of the robustness
#              simulation study,
#   Section 4: analysis of the simulation studies.
#
# The plotting helpers are in plot_functions.R. marginal_means() requires
# {emmeans}; the rest use base R.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Population parameters of the validation simulation study
# ─────────────────────────────────────────────────────────────────────────────

# True SD of the sensitivities that gives reliability reliab with k trials per
# condition (sampling variance as in usdt_reliability())
true_sd <- function(reliab, k, mu_d, mu_c = 0, sd_c = 0,
                    criterion_free = TRUE) {

  uniroot(function(sd_d) {
    d    <- mu_d + sd_d * qnorm(ppoints(200))
    crit <- if (sd_c > 0) mu_c + sd_c * qnorm(ppoints(20)) else mu_c
    grid <- expand.grid(d = d, crit = crit)

    eta_s <-  grid$d / 2 - grid$crit
    eta_n <- -grid$d / 2 - grid$crit
    w_s   <- k * dnorm(eta_s)^2 / (pnorm(eta_s) * (1 - pnorm(eta_s)))
    w_n   <- k * dnorm(eta_n)^2 / (pnorm(eta_n) * (1 - pnorm(eta_n)))

    error <- if (criterion_free) mean(1 / w_s + 1 / w_n) else mean(4 / (w_s + w_n))
    sd_d^2 / (sd_d^2 + error) - reliab
  }, interval = c(.001, 1.5), tol = 1e-8)$root
}

# True parameters of one condition (indirect SD and mean are iterated, since
# each depends on the other)
population_parameters <- function(reliab, trials, rho, intercept,
                                  mu_dD, mu_cD, sd_cD) {

  sd_dD <- true_sd(reliab, trials, mu_dD, mu_c = mu_cD, sd_c = sd_cD)

  sd_dI <- sd_dD
  repeat {
    mu_dI  <- intercept + rho * sd_dI / sd_dD * mu_dD
    sd_old <- sd_dI
    sd_dI  <- true_sd(reliab, trials, mu_dI, criterion_free = FALSE)
    if (abs(sd_dI - sd_old) < 1e-8) break
  }

  slope <- rho * sd_dI / sd_dD
  c(sd_dD  = sd_dD,                            # SD of true direct sensitivity
    sd_dI  = sd_dI,                            # SD of true indirect sensitivity
    mu_dI  = intercept + slope * mu_dD,        # mean true indirect sensitivity
    slope  = slope,                            # latent slope of d_I on d_D
    sd_res = sd_dI * sqrt(1 - rho^2),          # residual SD of d_I given d_D
    cov_d  = rho * sd_dD * sd_dI,              # covariance of the sensitivities
    diff   = intercept + slope * mu_dD - mu_dD) # group-level difference (H1)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Population parameters and reliability of the robustness study
# ─────────────────────────────────────────────────────────────────────────────

# Generating values and true hsdt parameters of one condition, for the
# Klauer et al. (1998) and Goldstein et al. (2022) models
robustness_parameters <- function(DGM, prop_unaware, rho, intercept,
                                  sd_dI, scale = 2) {

  if (DGM == "Klauer") {
    # d'_D = max(0, N), N ~ Normal(m_K, s_K), with P(d'_D = 0) = prop_unaware
    s_K <- scale
    m_K <- -s_K * qnorm(prop_unaware)
    a   <- m_K / s_K
    mu_dD  <- m_K * pnorm(a) + s_K * dnorm(a)
    var_dD <- (m_K^2 + s_K^2) * pnorm(a) + m_K * s_K * dnorm(a) - mu_dD^2

    # d'_I = intercept + beta * d'_D + e
    beta  <- rho * sd_dI / sqrt(var_dD)
    mu_x  <- mu_dD
    cov_d <- beta * var_dD
  } else {
    # theta = .5 (unaware) or U(.5, 1); d'_D = 2 * qnorm(theta) = 2|Z|
    m_K <- s_K <- NA
    p      <- 1 - prop_unaware
    mu_dD  <- p * 2 * sqrt(2 / pi)
    var_dD <- p * 4 - mu_dD^2

    # d'_I = intercept + beta * (theta - .5) + e, and E[2|Z|(Phi(|Z|) - .5)] = 1/sqrt(pi)
    mu_x  <- p / 4
    var_x <- p / 12 - mu_x^2
    beta  <- rho * sd_dI / sqrt(var_x)
    cov_d <- beta * (p / sqrt(pi) - mu_dD * mu_x)
  }

  mu_dI <- intercept + beta * mu_x
  slope <- cov_d / var_dD
  c(m_K            = m_K,                       # latent mean (Klauer)
    s_K            = s_K,                       # latent SD (Klauer)
    beta           = beta,                      # generating slope
    sd_e           = sd_dI * sqrt(1 - rho^2),   # residual SD of d'_I
    mu_dD          = mu_dD,                     # mean true direct sensitivity
    sd_dD          = sqrt(var_dD),              # SD of true direct sensitivity
    mu_dI          = mu_dI,                     # mean true indirect sensitivity
    cov_d          = cov_d,                     # covariance of the sensitivities
    cor_d          = cov_d / (sqrt(var_dD) * sd_dI), # correlation (H2)
    slope_hsdt     = slope,                     # slope of the linear projection
    intercept_hsdt = mu_dI - slope * mu_dD,     # hsdt intercept (H3)
    diff           = mu_dI - mu_dD)             # group-level difference (H1)
}

# Observed sensitivity with the log-linear correction (Hautus, 1995), with the
# criterion estimated (direct task) or fixed at zero (indirect task, where hits
# and correct rejections both have probability pnorm(d / 2))
observed_dprime <- function(hits, fas, k, criterion_free = TRUE) {
  if (criterion_free) {
    qnorm((hits + .5) / (k + 1)) - qnorm((fas + .5) / (k + 1))
  } else {
    2 * qnorm((hits + k - fas + .5) / (2 * k + 1))
  }
}

# Reliability of the observed sensitivities of one large sample, as the squared
# correlation between true and observed d'. The expected sampling variance used
# in the validation study (true_sd()) diverges for the direct task of the
# robustness study, because the variance of d' grows without bound for
# near-perfect participants. The indirect task is binary with criterion zero
# ("counts"), or normal RTs with unit SD split at each participant's median
# ("rt")
empirical_reliability <- function(d_D, d_I, c_D, k, indirect = c("counts", "rt")) {
  indirect <- match.arg(indirect)
  n  <- length(d_D)
  hD <- rbinom(n, k, pnorm( d_D / 2 - c_D))
  fD <- rbinom(n, k, pnorm(-d_D / 2 - c_D))
  if (indirect == "counts") {
    hI <- rbinom(n, k, pnorm( d_I / 2))
    fI <- rbinom(n, k, pnorm(-d_I / 2))
  } else {
    rt_s <- matrix(rnorm(n * k, -d_I / 2, 1), n)   # signal: faster by d_I
    rt_n <- matrix(rnorm(n * k,  d_I / 2, 1), n)
    med  <- apply(cbind(rt_s, rt_n), 1, median)
    hI   <- rowSums(rt_s < med)
    fI   <- rowSums(rt_n < med)
  }
  c(reliab_D = cor(d_D, observed_dprime(hD, fD, k))^2,
    reliab_I = cor(d_I, observed_dprime(hI, fI, k, criterion_free = FALSE))^2)
}

# Reliability of one condition of the robustness study, from n participants
# generated as in uSDT_robustness_sim.R (it does not depend on sample size)
robustness_reliability <- function(DGM, prop_unaware, rho, intercept, trials,
                                   sd_dI = .3, mu_cD = 0, sd_cD = .3, n = 1e5) {
  pars <- robustness_parameters(DGM, prop_unaware, rho, intercept, sd_dI)
  if (DGM == "Klauer") {
    d_D <- pmax(0, rnorm(n, pars[["m_K"]], pars[["s_K"]]))
    x   <- d_D
  } else {
    theta <- ifelse(runif(n) < prop_unaware, .5, runif(n, .5, 1))
    d_D   <- 2 * qnorm(theta)
    x     <- theta - .5
  }
  d_I <- intercept + pars[["beta"]] * x + rnorm(n, 0, pars[["sd_e"]])
  c_D <- rnorm(n, mu_cD, sd_cD)
  empirical_reliability(d_D, d_I, c_D, trials,
                        indirect = if (DGM == "Klauer") "counts" else "rt")
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Analysis of the simulation studies
# ─────────────────────────────────────────────────────────────────────────────

# Data of one parameter and measure, with the design factors as factors
data_split <- function(df_long, par, measure) {
  df <- as.data.frame(df_long[which(df_long$parameter == par &
                                      df_long$measure == measure), ])
  for (j in setdiff(colnames(df), "value")) df[, j] <- as.factor(df[, j])
  df
}

# Marginal means per design factor (bold: coverage outside [.925, .975] or
# relative bias of 10% or more)
marginal_means <- function(anova_fit, highlight_crit = "none", format = "latex") {
  factors <- attr(terms(anova_fit), "term.labels")
  factors <- factors[!grepl(":", factors)]
  table_data <- do.call(rbind, lapply(factors, function(f) {
    em <- as.data.frame(suppressMessages(emmeans::emmeans(anova_fit, specs = f)))
    data.frame(Factor_Name = f, Level = as.character(em[[f]]), value = em$emmean)
  }))
  flag <- switch(highlight_crit,
                 coverage = table_data$value < .925 | table_data$value > .975,
                 relative = abs(table_data$value) >= .10,
                 rep(FALSE, nrow(table_data)))
  txt  <- sprintf("%.3f", table_data$value)
  bold <- if (format == "latex") paste0("\\textbf{", txt, "}") else paste0("<b>", txt, "</b>")
  table_data$value <- ifelse(flag, bold, txt)
  table_data
}

# Lin's concordance correlation, with its precision (r) and accuracy (C_b)
ccc <- function(truth, estimate) {
  m_x <- mean(truth)
  m_y <- mean(estimate)
  s_x <- sqrt(mean((truth - m_x)^2))
  s_y <- sqrt(mean((estimate - m_y)^2))
  r   <- cor(truth, estimate)
  C_b <- 2 / (s_x / s_y + s_y / s_x + (m_x - m_y)^2 / (s_x * s_y))
  c(CCC = r * C_b, r = r, C_b = C_b)
}

# ─────────────────────────────────────────────────────────────────────────────
