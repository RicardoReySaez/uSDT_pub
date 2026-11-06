# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : uSDT_robustness_sim.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 28-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Simulation study on the robustness of the multilevel SDT model of uSDT when
# data are generated from two other population models:
#
#   Klauer et al. (1998): true direct sensitivity is a normal censored at
#     zero, and true indirect sensitivity is linear in it. Binary responses in
#     both tasks.
#   Goldstein et al. (2022): a share of participants is fully unaware, and the
#     rest have accuracies from U(.5, 1). True indirect sensitivity is linear
#     in accuracy, not in d'. Indirect task measured with response times,
#     split at the median by uSDT.
#
# Both models get residual variance in the indirect sensitivity and criterion
# variability in the direct task. The intercept (H3) is compared with two
# true values: the intercept of each model (alpha) and the hsdt intercept of
# the population (they differ in Goldstein's model when the slope is not 0).

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────
# Libraries necessary for the script to function
library(SimDesign)
library(uSDT)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# robustness_parameters() (Section 3)
source("R scripts/R functions/simulation_functions.R")

# Fixed population parameters
mu_cD <- 0   # mean direct criterion
sd_cD <- .3  # SD of direct criterion (indirect criterion fixed at zero)
sd_dI <- .3  # SD of true indirect sensitivity

# Design conditions
Design <- createDesign(
  DGM          = c("Klauer", "Goldstein"),
  N            = c(30, 100, 200),
  trials       = c(10, 30, 50),
  prop_unaware = c(.2, .5, .8),
  rho          = c(0, .2, .4, .6),
  intercept    = c(0, .15, .30))

# Generating values and true population parameters of each condition
true_pars <- t(mapply(robustness_parameters,
                      DGM          = Design$DGM,
                      prop_unaware = Design$prop_unaware,
                      rho          = Design$rho,
                      intercept    = Design$intercept,
                      MoreArgs     = list(sd_dI = sd_dI)))
Design$mu_cD <- mu_cD
Design$sd_cD <- sd_cD
Design$sd_dI <- sd_dI
for (p in colnames(true_pars)) Design[[p]] <- true_pars[, p]

# Folder of the results
dir.create("Results/Simulation studies/Raw/Robustness study", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: data-generating function
# ─────────────────────────────────────────────────────────────────────────────

# condition <- Design[1, ]
Generate <- function(condition, fixed_objects) {
  Attach(condition)
  # True direct sensitivities, and the predictor of the indirect ones: d'
  # (Klauer) or accuracy above chance (Goldstein)
  if (DGM == "Klauer") {
    d_D <- pmax(0, rnorm(N, m_K, s_K))
    x   <- d_D
  } else {
    theta <- ifelse(runif(N) < prop_unaware, .5, runif(N, .5, 1))
    d_D   <- 2 * qnorm(theta)
    x     <- theta - .5
  }

  # True indirect sensitivities
  d_I <- intercept + beta * x + rnorm(N, 0, sd_e)

  # Direct criterion (indirect criterion fixed at zero)
  c_D <- rnorm(N, mu_cD, sd_cD)

  # Direct task: hit and false-alarm counts, in the format of usdt_data_tasks
  direct <- data.frame(
    subj      = rep(1:N, 2),
    condition = rep(c("signal", "noise"), each = N),
    yes       = c(rbinom(N, trials, pnorm( d_D / 2 - c_D)),
                  rbinom(N, trials, pnorm(-d_D / 2 - c_D))),
    n         = trials)

  # Indirect task: counts with criterion zero (Klauer), or standardized RTs,
  # faster by d_I in the signal condition (Goldstein)
  if (DGM == "Klauer") {
    indirect <- data.frame(
      subj      = rep(1:N, 2),
      condition = rep(c("signal", "noise"), each = N),
      yes       = c(rbinom(N, trials, pnorm( d_I / 2)),
                    rbinom(N, trials, pnorm(-d_I / 2))),
      n         = trials)
  } else {
    cond     <- rep(rep(c(.5, -.5), each = trials), N)
    indirect <- data.frame(
      subj      = rep(1:N, each = 2 * trials),
      condition = ifelse(cond > 0, "signal", "noise"),
      rt        = rnorm(2 * N * trials, -rep(d_I, each = 2 * trials) * cond, 1))
  }

  dat <- list(direct = direct, indirect = indirect)
  dat
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: data-analysis function
# ─────────────────────────────────────────────────────────────────────────────

Analyse <- function(condition, dat, fixed_objects) {
  # Data in uSDT format
  if (condition$DGM == "Klauer") {
    d <- uSDT::usdt_data_tasks(
      direct           = dat$direct,
      indirect         = dat$indirect,
      subject_col      = "subj",
      condition_col    = "condition",
      condition_levels = c(signal = "signal", noise = "noise"),
      successes_col    = "yes",
      trials_col       = "n",
      successes_type   = "counts")

    # TEMPORARY: hsdt() only fixes the indirect criterion at zero when uSDT did
    # Meyen's split itself, so the flag is set by hand
    d$meta$criterion_zero[["indirect"]] <- TRUE
  } else {
    # Indirect RTs split at each participant's median by uSDT
    d <- uSDT::usdt_data_tasks(
      direct           = dat$direct,
      indirect         = dat$indirect,
      subject_col      = "subj",
      condition_col    = "condition",
      condition_levels = c(signal = "signal", noise = "noise"),
      successes_col    = list(direct = "yes", indirect = NULL),
      trials_col       = list(direct = "n",   indirect = NULL),
      successes_type   = "counts",
      response_col     = list(direct = NULL, indirect = "rt"),
      response_levels  = list(direct = NULL, indirect = c(signal = "faster", noise = "slower")),
      dichotomize      = list(direct = FALSE, indirect = TRUE))
  }

  # Estimate multilevel SDT model
  time_fit <- system.time({
    fit <- tryCatch(suppressWarnings(suppressMessages(uSDT::hsdt(d))),
                    error = function(e) NULL)
  })

  # Estimate, SE, CI and p value of each parameter. A failed fit keeps NA, so
  # that it is counted instead of replaced by new data
  pars <- c("mu_dD", "mu_dI", "mu_cD", "sd_dD", "sd_dI", "sd_cD", "cov_d",
            "diff", "rho", "intercept", "slope")
  out  <- matrix(NA_real_, length(pars), 5,
                 dimnames = list(pars, c("est", "se", "lo", "hi", "p")))

  # Diagnostics. singular is lme4's boundary (singular) fit warning; a
  # hypothesis is testable only if it has SE and CI
  diagnostics <- c(failed = 1, converged = NA, singular = NA,
                   testable_H1 = 0, testable_H2 = 0, testable_H3 = 0)

  if (!is.null(fit)) {
    z <- qnorm(.975)

    # Fixed effects (the criterion is fitted with the opposite sign)
    out[c("mu_cD", "mu_dD", "mu_dI"), "est"] <-
      lme4::fixef(fit$fit)[c("c_D", "d_D", "d_I")] * c(-1, 1, 1)
    out[c("mu_cD", "mu_dD", "mu_dI"), "se"]  <-
      sqrt(diag(as.matrix(vcov(fit$fit))))[c("c_D", "d_D", "d_I")]

    # Random-effect SDs
    re <- fit$pars$re
    out[c("sd_cD", "sd_dD", "sd_dI"), "est"] <- re$sd[match(c("c_D", "d_D", "d_I"), re$term)]
    out[c("sd_cD", "sd_dD", "sd_dI"), "se"]  <- re$se[match(c("c_D", "d_D", "d_I"), re$term)]

    # Covariance between sensitivities
    out["cov_d", "est"] <- fit$pars$est[["s_DI"]]
    out["cov_d", "se"]  <- suppressWarnings(sqrt(fit$pars$vcov["s_DI", "s_DI"]))

    # Wald CIs and p values (log scale for SDs, as in summary.hsdt)
    wald <- c("mu_cD", "mu_dD", "mu_dI", "cov_d")
    sds  <- c("sd_cD", "sd_dD", "sd_dI")
    out[wald, "lo"] <- out[wald, "est"] - z * out[wald, "se"]
    out[wald, "hi"] <- out[wald, "est"] + z * out[wald, "se"]
    out[wald, "p"]  <- 2 * pnorm(-abs(out[wald, "est"] / out[wald, "se"]))
    out[sds, "lo"]  <- exp(log(out[sds, "est"]) - z * out[sds, "se"] / out[sds, "est"])
    out[sds, "hi"]  <- exp(log(out[sds, "est"]) + z * out[sds, "se"] / out[sds, "est"])

    # Three hypotheses
    tests <- fit$tests[match(c("d'(indirect) - d'(direct)", "correlation",
                               "intercept", "slope"), fit$tests$term), ]
    out[c("diff", "rho", "intercept", "slope"), ] <-
      as.matrix(tests[, c("estimate", "se", "conf.low", "conf.high", "p.value")])

    diagnostics <- c(failed      = 0,
                     converged   = fit$diagnostics$converged,
                     singular    = fit$diagnostics$singular,
                     testable_H1 = tests$status[1] %in% "ok",
                     testable_H2 = tests$status[2] %in% "ok",
                     testable_H3 = tests$status[3] %in% "ok")
  }

  # Return results: <parameter>_<statistic>, diagnostics and time
  res <- c(setNames(c(t(out)), paste(rep(pars, each = 5), colnames(out), sep = "_")),
           diagnostics,
           time_fit = time_fit[["elapsed"]])

  # Return results
  return(res)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: summarise function
# ─────────────────────────────────────────────────────────────────────────────

Summarise <- function(condition, results, fixed_objects) {
  Attach(condition)
  # True values. The intercept estimate is compared with the hsdt intercept
  # of the population (intercept) and with the model's own (alpha)
  truth <- c(mu_dD = mu_dD, mu_dI = mu_dI, mu_cD = mu_cD, sd_dD = sd_dD,
             sd_dI = sd_dI, sd_cD = sd_cD, cov_d = cov_d, diff = diff,
             rho = cor_d, intercept = intercept_hsdt, slope = slope_hsdt,
             alpha = intercept)

  # Bias, relative bias and RMSE over the replications with an estimate;
  # coverage and MSRSE over those with CI and SE
  recovery <- unlist(lapply(names(truth), function(p) {
    q   <- if (p == "alpha") "intercept" else p
    est <- results[[paste0(q, "_est")]]
    se  <- results[[paste0(q, "_se")]]
    CIs <- cbind(results[[paste0(q, "_lo")]], results[[paste0(q, "_hi")]])
    ok    <- !is.na(est)
    se_ok <- ok & !is.na(se)
    ci_ok <- complete.cases(CIs)
    res <- c(bias        = if (any(ok)) bias(est[ok], truth[[p]]) else NA,
             rel_bias    = if (any(ok) && truth[[p]] != 0) bias(est[ok], truth[[p]], type = "relative") else NA,
             RMSE        = if (any(ok)) RMSE(est[ok], truth[[p]]) else NA,
             coverage    = if (any(ci_ok)) ECR(CIs[ci_ok, , drop = FALSE], truth[[p]]) else NA,
             MSRSE       = if (sum(se_ok) > 1) MSRSE(cbind(se[se_ok]), cbind(est[se_ok])) else NA,
             missing_est = mean(!ok),
             missing_ci  = mean(!ci_ok))
    setNames(res, paste0(p, "_", names(res)))
  }))

  # Rejection rates of the three hypotheses over the testable replications
  # (type I error when the truth is zero, power otherwise)
  rejection <- sapply(c("diff", "rho", "intercept", "slope"), function(p) {
    pv <- results[[paste0(p, "_p")]]
    if (any(!is.na(pv))) EDR(pv[!is.na(pv)], alpha = .05) else NA
  })
  names(rejection) <- paste0(names(rejection), "_rejection")

  # Boundary: singular fits, hypotheses without SE (not testable), and which
  # variance is estimated at zero or correlation at +-1
  boundary <- c(
    singular      = mean(results$singular, na.rm = TRUE),
    untestable_H1 = 1 - mean(results$testable_H1),
    untestable_H2 = 1 - mean(results$testable_H2),
    untestable_H3 = 1 - mean(results$testable_H3),
    zero_sd_dD    = mean(!is.na(results$sd_dD_est) & results$sd_dD_est < 1e-4),
    zero_sd_dI    = mean(!is.na(results$sd_dI_est) & results$sd_dI_est < 1e-4),
    zero_sd_cD    = mean(!is.na(results$sd_cD_est) & results$sd_cD_est < 1e-4),
    one_rho       = mean(!is.na(results$rho_est) & abs(results$rho_est) > 1 - 1e-6))

  # Convergence
  fits <- c(failed    = mean(results$failed),
            converged = mean(results$converged, na.rm = TRUE),
            time_fit  = mean(results$time_fit))

  c(recovery, rejection, boundary, fits)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: simulation
# ─────────────────────────────────────────────────────────────────────────────

# Global seeds per condition
global_seeds <- genSeeds(design = Design, iseed = 2026)

# allow_na and allow_nan keep boundary and failed fits (NA or NaN estimates)
# instead of redrawing data
res <- runSimulation(
  design         = Design,
  replications   = 1000,
  generate       = Generate,
  analyse        = Analyse,
  summarise      = Summarise,
  parallel       = TRUE,
  store_results  = TRUE,
  save           = TRUE,
  max_errors     = 1000L,
  ncores         = 12,
  seed           = global_seeds,
  save_details   = list(tmpfilename = "tmp_uSDT_robustness_simres.qs2",
                        out_rootdir = "Results/Simulation studies/Raw/Robustness study"),
  control        = list(store_Random.seeds = TRUE, allow_na = TRUE,
                        allow_nan = TRUE),
  packages       = c("uSDT", "lme4"),
  filename       = "uSDT_robustness_simres.rds"
)

# Store again results as rdata
save(res, file = "Results/Simulation studies/Raw/Robustness study/uSDT_robustness_simres.rdata")

# ─────────────────────────────────────────────────────────────────────────────


