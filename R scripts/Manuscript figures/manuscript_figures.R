# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : manuscript_figures.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 02-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Make Figures 1 to 4 of the main text in Figures/Manuscript/. Self-contained:
# it needs only ggplot2, patchwork and uSDT.
# Figures 3 and 4 share one population, awareness-dependent by construction
# (zero latent intercept); observed scores come from simulated trials.
# Run from the folder containing uSDT_pub.Rproj.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Packages and settings
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
for (pkg in c("uSDT", "patchwork")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("The ", pkg, " package is needed.", call. = FALSE)
  }
}
if (!file.exists("uSDT_pub.Rproj")) {
  stop("Run this script from the folder holding uSDT_pub.Rproj.", call. = FALSE)
}
dir.create("Figures/Manuscript", recursive = TRUE, showWarnings = FALSE)

# Seed of the simulations, and of the sample of 40 participants shown in
# Figures 3 and 4
seed        <- 2025
sample_seed <- 9637

# True sensitivities: d_D ~ Normal(mu_d, sd_d) and d_I = slope * d_D + error,
# so the latent intercept is zero and the mean difference is −.15
pars <- within(list(mu_d = 0.50, sd_d = 0.45, slope = 0.70, sd_resid = 0.25), {
  mu_i      <- slope * mu_d
  sd_i      <- sqrt(slope^2 * sd_d^2 + sd_resid^2)
  rho       <- slope * sd_d / sd_i
  intercept <- mu_i - slope * mu_d
})

# Trials per task, half signal and half noise. The direct criterion varies
# across participants; the indirect one is zero, as after Meyen et al.'s
# median split.
designs <- list(moderate = c(direct = 80, indirect = 80),
                few      = c(direct = 20, indirect = 20))
criterion_sd <- 0.30

n_sample <- 40   # participants in the sample
n_limit  <- 2e5  # stands in for the population

# Figure 1: power of the direct test for a small per-trial effect
ita <- list(n            = 30,
            d            = 0.15,
            sd_between   = 0.10,
            trials_grid  = unique(c(1, round(2^seq(1, 9, by = 0.25) / 2) * 2)),
            reps         = 10000,
            reps_typical = 20000)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Functions
# ─────────────────────────────────────────────────────────────────────────────

# True sensitivities of n participants
simulate_true <- function(n) {
  d_direct <- rnorm(n, pars$mu_d, pars$sd_d)
  list(d_direct   = d_direct,
       d_indirect = pars$slope * d_direct + rnorm(n, 0, pars$sd_resid))
}

# Observed d' from simulated trials: equal-variance SDT with a criterion drawn
# from Normal(0, criterion_sd), and the log-linear correction (Hautus, 1995).
# The counts are kept for the multilevel model.
observe <- function(true_d, n_trials, criterion_sd) {
  k    <- n_trials / 2
  crit <- rnorm(length(true_d), 0, criterion_sd)
  hits <- rbinom(length(true_d), k, pnorm( true_d / 2 - crit))
  fas  <- rbinom(length(true_d), k, pnorm(-true_d / 2 - crit))
  d    <- qnorm((hits + 0.5) / (k + 1)) - qnorm((fas + 0.5) / (k + 1))
  structure(d, hits = hits, fas = fas, k = k)
}

observe_both <- function(true, trials) {
  list(d_direct   = observe(true$d_direct,   trials[["direct"]],   criterion_sd),
       d_indirect = observe(true$d_indirect, trials[["indirect"]], 0))
}

# Means, SDs, correlation and regression of indirect on direct sensitivity
moments <- function(d_direct, d_indirect) {
  slope <- cov(d_direct, d_indirect) / var(d_direct)
  list(mu_d      = mean(d_direct),
       mu_i      = mean(d_indirect),
       sd_d      = sd(d_direct),
       sd_i      = sd(d_indirect),
       rho       = cor(d_direct, d_indirect),
       slope     = slope,
       intercept = mean(d_indirect) - slope * mean(d_direct))
}

# Two-step estimates in one sample: regression intercept (95% CI) and slope,
# and correlation (95% CI from Fisher's z)
two_step <- function(d_direct, d_indirect) {
  fit <- lm(d_indirect ~ d_direct)
  ci  <- confint(fit)
  r   <- cor(d_direct, d_indirect)
  se  <- 1 / sqrt(length(d_direct) - 3)
  c(intercept    = coef(fit)[[1]],
    intercept_lo = ci[1, 1],
    intercept_hi = ci[1, 2],
    slope        = coef(fit)[[2]],
    r            = r,
    r_lo         = tanh(atanh(r) - qnorm(0.975) * se),
    r_hi         = tanh(atanh(r) + qnorm(0.975) * se))
}

# Multilevel SDT model (uSDT::hsdt) fitted to the trials behind a sample's
# observed d'. The indirect criterion is zero by construction, so it is fixed
# at zero by hand. Returns the latent regression and correlation with 95%
# CIs, and each participant's conditional sensitivities (fixed effect plus
# conditional mode).
fit_multilevel <- function(obs) {
  counts <- function(d) {
    data.frame(subj      = rep(seq_along(d), 2),
               condition = rep(c("signal", "noise"), each = length(d)),
               yes       = c(attr(d, "hits"), attr(d, "fas")),
               n         = attr(d, "k"))
  }
  data <- uSDT::usdt_data_tasks(
    awareness        = counts(obs$d_direct),
    influence        = counts(obs$d_indirect),
    subject_col      = "subj",
    condition_col    = "condition",
    condition_levels = c(signal = "signal", noise = "noise"),
    successes_col    = "yes",
    trials_col       = "n",
    successes_type   = "counts"
  )
  data$meta$criterion_zero[["influence"]] <- TRUE

  m     <- suppressWarnings(suppressMessages(uSDT::hsdt(data)))
  tests <- m$tests
  rownames(tests) <- tests$term
  re <- lme4::ranef(m$fit)$subj[as.character(seq_along(obs$d_direct)), ]
  fe <- lme4::fixef(m$fit)

  list(estimate    = c(intercept    = tests["intercept",   "estimate"],
                       intercept_lo = tests["intercept",   "conf.low"],
                       intercept_hi = tests["intercept",   "conf.high"],
                       r            = tests["correlation", "estimate"],
                       r_lo         = tests["correlation", "conf.low"],
                       r_hi         = tests["correlation", "conf.high"],
                       slope        = tests["slope",       "estimate"]),
       conditional = data.frame(d_direct   = fe[["d_A"]] + re$d_A,
                                d_indirect = fe[["d_I"]] + re$d_I))
}

# One named value from each element of a list
pick <- function(x, name) unname(sapply(x, `[[`, name))

# Number formats: true minus sign; correlations without a leading zero
fmt_num <- function(x) sub("^-", "−", sprintf("%.2f", x))
fmt_r   <- function(x) sub("^-", "−", sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x)))
ci_text <- function(est, name, fmt) {
  sprintf("%s [%s, %s]", fmt(est[name]), fmt(est[paste0(name, "_lo")]),
          fmt(est[paste0(name, "_hi")]))
}

# Shared look and axis titles
theme_fig <- theme_classic(base_size = 9) +
  theme(strip.background = element_blank(),
        strip.text       = element_text(hjust = 0, face = "bold", size = 8.5),
        panel.border     = element_rect(fill = NA, colour = "grey80"),
        axis.title       = element_text(size = 8.5))
labs_d <- labs(x = quote("Awareness sensitivity, " * italic(d) * minute[A]),
               y = quote("Influence sensitivity, " * italic(d) * minute[I]))
zero_line_colour <- "steelblue"

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Simulate the data of Figures 3 and 4
# ─────────────────────────────────────────────────────────────────────────────

# Population: true values, and observed-score moments under each trial count
# from one large sample, observed with separate seeds
set.seed(seed)
big <- simulate_true(n_limit)

observe_population <- function(trials, s) {
  set.seed(s)
  obs <- observe_both(big, trials)
  moments(obs$d_direct, obs$d_indirect)
}
population <- list(true     = pars,
                   moderate = observe_population(designs$moderate, seed + 2),
                   few      = observe_population(designs$few,      seed + 1))

# Sample: the same participants under both trial counts, analysed in two
# steps (Figure 3) and with the multilevel model (Figure 4)
set.seed(sample_seed)
latent       <- simulate_true(n_sample)
observed     <- lapply(designs[c("moderate", "few")], observe_both, true = latent)
samples      <- c(list(true = latent), lapply(observed, lapply, as.vector))
two_step_est <- lapply(samples, function(s) two_step(s$d_direct, s$d_indirect))
multilevel   <- lapply(observed, fit_multilevel)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Figure 1, power of one sensitivity as a function of trials
# ─────────────────────────────────────────────────────────────────────────────

# Each participant classifies k trials, correct with probability pnorm(d / 2)
# and d ~ Normal(ita$d, ita$sd_between); the proportions correct are tested
# against chance (.5) with a one-sample t test, as in the original analyses.
power_direct <- function(k, reps) {
  n   <- ita$n
  d_i <- matrix(rnorm(n * reps, ita$d, ita$sd_between), n, reps)
  acc <- matrix(rbinom(n * reps, k, pnorm(d_i / 2)), n, reps) / k - 0.5
  t   <- colMeans(acc) / (apply(acc, 2, sd) / sqrt(n))
  p   <- 2 * pt(-abs(t), df = n - 1)
  mean(!is.na(p) & p < .05)
}

set.seed(seed)
power_curve <- data.frame(k     = ita$trials_grid,
                          power = sapply(ita$trials_grid, power_direct,
                                         reps = ita$reps))

# Typical awareness (5, 20 trials) and influence (200 trials) designs, with
# their labels kept away from the curve
typical <- data.frame(k       = c(5, 20, 200),
                      nudge_x = c(0.85, 0.85, 1.12),
                      nudge_y = c(0.08, 0.08, -0.08),
                      hjust   = c(1, 1, 0))
typical$power <- sapply(typical$k, power_direct, reps = ita$reps_typical)
typical$label <- sprintf("%d trials: %.0f%%", typical$k, 100 * typical$power)

ranges <- data.frame(xmin  = c(1, 100), xmax = c(20, 400),
                     label = c("Typical awareness tasks",
                               "Typical influence tasks"))

fig1 <- ggplot(power_curve, aes(k, power)) +
  geom_rect(data = ranges, aes(xmin = xmin, xmax = xmax, ymin = -Inf,
                               ymax = Inf),
            fill = "grey93", inherit.aes = FALSE) +
  geom_text(data = ranges, aes(x = sqrt(xmin * xmax), y = 1.06,
                               label = label),
            size = 2.5, fontface = "italic", inherit.aes = FALSE) +
  geom_hline(yintercept = 0.80, linetype = "13", colour = "grey50",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(data = typical, shape = 21, size = 1.7, fill = "white",
             stroke = 0.7) +
  geom_text(data = typical,
            aes(x = k * nudge_x, y = power + nudge_y, label = label,
                hjust = hjust),
            size = 2.5) +
  scale_x_continuous(trans = "log2", breaks = c(1, 4, 16, 64, 256),
                     expand = expansion(mult = c(0.02, 0.09))) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1.1), clip = "off") +
  labs(x = "Number of trials per participant",
       y = "Probability of a significant result") +
  theme_fig

ggsave("Figures/Manuscript/Figure_1.png", fig1, width = 4.8, height = 3.1)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Figure 2, one RT effect and two sensitivities
# ─────────────────────────────────────────────────────────────────────────────

# Same 50-ms effect, different trial-to-trial SDs. Top: RT distributions in
# ms. Bottom: standardized distributions, whose means differ by d'. These are
# population normals, not simulated samples.
congruent <- 600
effect    <- 50
rt_sd     <- c(25, 125)
dprime    <- effect / rt_sd

panels_rt <- c("A. Participant 1, in milliseconds",
               "B. Participant 2, in milliseconds",
               "C. Participant 1, in trial-to-trial SDs",
               "D. Participant 2, in trial-to-trial SDs")
conditions_rt <- c("Congruent trials", "Incongruent trials")

rt_density <- function(p, scale) {
  if (scale == "ms") {
    x <- seq(200, 1050, length.out = 801)
    means <- congruent + c(0, effect)
    s <- rt_sd[p]
  } else {
    x <- seq(-4, 6, length.out = 801)
    means <- c(0, dprime[p])
    s <- 1
  }
  panel <- factor(panels_rt[p + if (scale == "ms") 0 else 2], levels = panels_rt)
  do.call(rbind, lapply(1:2, function(k)
    data.frame(panel     = panel,
               condition = factor(conditions_rt[k], levels = conditions_rt),
               x         = x,
               density   = dnorm(x, means[k], s))))
}

density_rt <- do.call(rbind, c(lapply(1:2, rt_density, scale = "ms"),
                               lapply(1:2, rt_density, scale = "sd")))

labels_rt <- data.frame(
  panel = factor(panels_rt, levels = panels_rt),
  label = c(sprintf("Effect = %d ms\nSD = %d ms", effect, rt_sd),
            sprintf("d′ = %s", fmt_num(dprime)))
)

# Separate rows allow different x-axis titles.
row_rt <- function(rows, x_title) {
  keep <- panels_rt[rows]
  ggplot(density_rt[density_rt$panel %in% keep, ], aes(x, density)) +
    geom_area(aes(fill = condition), position = "identity", alpha = 0.5) +
    geom_line(aes(linetype = condition), linewidth = 0.6) +
    geom_text(data = labels_rt[labels_rt$panel %in% keep, ],
              aes(x = Inf, y = Inf, label = label),
              hjust = 1.08, vjust = 1.3, size = 2.6, lineheight = 0.95,
              inherit.aes = FALSE) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("grey60", "grey90"), name = NULL) +
    scale_linetype_manual(values = c("solid", "22"), name = NULL) +
    scale_x_continuous(labels = function(x) sub("^-", "−", x)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(x = x_title, y = "Density") +
    theme_fig +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank())
}

fig2 <- patchwork::wrap_plots(row_rt(1:2, "Response time (ms)"),
                            row_rt(3:4, "z-score"), ncol = 1) +
  patchwork::plot_layout(guides = "collect", axis_titles = "collect_y") &
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.9, "cm"),
        legend.text      = element_text(size = 8))

ggsave("Figures/Manuscript/Figure_2.png", fig2, width = 6.0, height = 4.4)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Figure 3, participants by number of trials
# ─────────────────────────────────────────────────────────────────────────────

# Columns: infinite, moderate or few trials. Rows: population or 40 people.
# Individual sensitivities are computed without a multilevel model.
trials_label <- function(trials) sprintf("%d trials per task", trials[["direct"]])
sample_label <- sprintf("%d participants", n_sample)

panels_f2 <- c(
  "A. Infinite participants, infinite trials",
  paste0("B. Infinite participants, ", trials_label(designs$moderate)),
  paste0("C. Infinite participants, ", trials_label(designs$few)),
  paste0("D. ", sample_label, ", infinite trials"),
  paste0("E. ", sample_label, ", ", trials_label(designs$moderate)),
  paste0("F. ", sample_label, ", ", trials_label(designs$few))
)
as_panel    <- function(i) factor(panels_f2[i], levels = panels_f2)
axis_limits <- c(-2.2, 2.5)

# Top row: bivariate normal density inside its 95% region, faded to white at
# the edge; sqrt keeps faint densities visible.
edge <- sqrt(exp(-qchisq(0.95, df = 2) / 2))
density_f2 <- do.call(rbind, lapply(1:3, function(i) {
  p     <- population[[i]]
  axis  <- seq(axis_limits[1], axis_limits[2], length.out = 400)
  grid  <- expand.grid(d_direct = axis, d_indirect = axis)
  cov   <- p$rho * p$sd_d * p$sd_i
  dist  <- mahalanobis(grid, c(p$mu_d, p$mu_i),
                       matrix(c(p$sd_d^2, cov, cov, p$sd_i^2), 2, 2))
  keep  <- dist <= qchisq(0.95, df = 2)
  data.frame(grid[keep, ],
             shade = (sqrt(exp(-dist[keep] / 2)) - edge) / (1 - edge),
             panel = as_panel(i))
}))

# Bottom row: the same participants at each trial count
points_f2 <- do.call(rbind, lapply(1:3, function(i)
  data.frame(samples[[i]], panel = as_panel(i + 3))))

# Regression line and intercept of each panel: population values on top,
# two-step estimates with 95% CIs below
est_f2   <- c(population, two_step_est)
lines_f2 <- data.frame(d_direct   = 0,
                       d_indirect = pick(est_f2, "intercept"),
                       intercept  = pick(est_f2, "intercept"),
                       slope      = pick(est_f2, "slope"),
                       lo         = c(NA, NA, NA, pick(two_step_est, "intercept_lo")),
                       hi         = c(NA, NA, NA, pick(two_step_est, "intercept_hi")),
                       panel      = as_panel(1:6))

# Solid in A: true line. Dashed elsewhere: true-line reference.
true_line <- data.frame(intercept = pars$intercept, slope = pars$slope,
                        panel = as_panel(1:6),
                        kind  = c("solid", rep("22", 5)),
                        ink   = c("black", rep("grey30", 5)),
                        width = c(0.7, rep(0.5, 5)))

labels_f2 <- data.frame(
  panel = as_panel(1:6),
  txt   = c(sprintf("Intercept = %s\nρ = %s",
                    fmt_num(pick(population, "intercept")),
                    fmt_r(pick(population, "rho"))),
            sprintf("Intercept = %s\nr = %s",
                    sapply(two_step_est, ci_text, "intercept", fmt_num),
                    sapply(two_step_est, ci_text, "r", fmt_r)))
)

fig3 <- ggplot(mapping = aes(d_direct, d_indirect)) +
  geom_raster(data = density_f2, aes(fill = shade), interpolate = TRUE) +
  scale_fill_gradient(low = "white", high = "grey40", limits = c(0, 1),
                      guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = zero_line_colour) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = zero_line_colour) +
  geom_point(data = points_f2, shape = 16, size = 1, colour = "grey45",
             alpha = 0.55) +
  geom_abline(data = true_line,
              aes(intercept = intercept, slope = slope, linetype = kind,
                  colour = ink, linewidth = width)) +
  scale_colour_identity() + scale_linewidth_identity() +
  scale_linetype_identity() +
  geom_abline(data = lines_f2[-1, ], aes(intercept = intercept, slope = slope),
              linewidth = 0.7) +
  geom_errorbar(data = lines_f2, aes(ymin = lo, ymax = hi),
                width = 0.18, linewidth = 0.55, na.rm = TRUE) +
  geom_point(data = lines_f2, shape = 21, size = 1.7,
             fill = "white", stroke = 0.7) +
  geom_text(data = labels_f2, aes(label = txt),
            x = axis_limits[2] - 0.05, y = axis_limits[2] - 0.05,
            hjust = 1, vjust = 1, size = 2.4, lineheight = 1.05) +
  facet_wrap(~panel, ncol = 3) +
  coord_fixed(xlim = axis_limits, ylim = axis_limits, expand = FALSE) +
  labs_d +
  theme_fig +
  theme(strip.text = element_text(size = 8))

ggsave("Figures/Manuscript/Figure_3.png", fig3, width = 7.6, height = 5.6)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 8: Figure 4, the multilevel model fitted to the same trials
# ─────────────────────────────────────────────────────────────────────────────

# Segments show shrinkage from observed to conditional sensitivities;
# solid = fitted latent line, dashed = true line.
# Our observed scores use the log-linear correction for everyone; uSDT's
# regression plot uses it only for extreme rates.
panels_f3 <- c(
  moderate = paste0("A. ", sample_label, ", ", trials_label(designs$moderate)),
  few      = paste0("B. ", sample_label, ", ", trials_label(designs$few))
)
as_panel_f3 <- function(nm) factor(unname(panels_f3[nm]), levels = panels_f3)

shrinkage_f3 <- do.call(rbind, lapply(names(multilevel), function(nm)
  data.frame(observed    = samples[[nm]],
             conditional = multilevel[[nm]]$conditional,
             panel       = as_panel_f3(nm))))

# Latent regression and intercept with 95% CI
estimates_f3 <- lapply(multilevel, `[[`, "estimate")
lines_f3 <- data.frame(d_direct   = 0,
                       d_indirect = pick(estimates_f3, "intercept"),
                       intercept  = pick(estimates_f3, "intercept"),
                       slope      = pick(estimates_f3, "slope"),
                       lo         = pick(estimates_f3, "intercept_lo"),
                       hi         = pick(estimates_f3, "intercept_hi"),
                       panel      = as_panel_f3(names(multilevel)))

# Hatted rho for model estimates, to distinguish it from observed r
labels_f3 <- data.frame(
  panel     = as_panel_f3(names(multilevel)),
  intercept = sprintf("Intercept = %s",
                      sapply(estimates_f3, ci_text, "intercept", fmt_num)),
  rho       = sprintf('hat("ρ") == "%s"',
                      sapply(estimates_f3, ci_text, "r", fmt_r))
)

fig4 <- ggplot(mapping = aes(d_direct, d_indirect)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = zero_line_colour) +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = zero_line_colour) +
  geom_segment(data = shrinkage_f3,
               aes(x = observed.d_direct, y = observed.d_indirect,
                   xend = conditional.d_direct, yend = conditional.d_indirect),
               colour = "grey78", linewidth = 0.25, inherit.aes = FALSE) +
  geom_point(data = shrinkage_f3,
             aes(x = observed.d_direct, y = observed.d_indirect),
             shape = 16, size = 1, colour = "grey45", alpha = 0.55,
             inherit.aes = FALSE) +
  geom_point(data = shrinkage_f3,
             aes(x = conditional.d_direct, y = conditional.d_indirect),
             shape = 16, size = 1, colour = "black", inherit.aes = FALSE) +
  geom_abline(intercept = pars$intercept, slope = pars$slope,
              linetype = "22", colour = "grey30", linewidth = 0.5) +
  geom_abline(data = lines_f3, aes(intercept = intercept, slope = slope),
              linewidth = 0.7) +
  geom_errorbar(data = lines_f3, aes(ymin = lo, ymax = hi),
                width = 0.18, linewidth = 0.55) +
  geom_point(data = lines_f3, shape = 21, size = 1.7,
             fill = "white", stroke = 0.7) +
  geom_text(data = labels_f3, aes(label = intercept),
            x = axis_limits[2] - 0.05, y = axis_limits[2] - 0.05,
            hjust = 1, vjust = 1, size = 2.4) +
  geom_text(data = labels_f3, aes(label = rho), parse = TRUE,
            x = axis_limits[2] - 0.05, y = axis_limits[2] - 0.05 - 0.18,
            hjust = 1, vjust = 1, size = 2.4) +
  facet_wrap(~panel, ncol = 2) +
  coord_fixed(xlim = axis_limits, ylim = axis_limits, expand = FALSE) +
  labs_d +
  theme_fig +
  theme(strip.text = element_text(size = 8))

ggsave("Figures/Manuscript/Figure_4.png", fig4, width = 6.0, height = 3.4)

# ─────────────────────────────────────────────────────────────────────────────
