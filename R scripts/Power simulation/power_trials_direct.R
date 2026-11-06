# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : power_trials_direct.R                                      ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 16-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Power for H1 (indirect minus direct sensitivity) across direct-task trials.
# Output: Figures/Manuscript/Figure_5.png
# Uses the multilevel SDT Wald test through uSDT::usdt_power().
# The approximation differed by at most .04 from 1,000 simulated studies
# per design. Use method = "simulation" to rerun that slower check.
# Run from the project root (uSDT_pub.Rproj).
# Last run: R 4.5.0; uSDT 0.1.0.9000 (dev, commit 795096f).

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

library(ggplot2)
if (!requireNamespace("uSDT", quietly = TRUE) ||
    !"usdt_power" %in% getNamespaceExports("uSDT")) {
  stop("This script needs a version of uSDT with usdt_power(): ",
       "remotes::install_github(\"RicardoReySaez/uSDT\", ref = \"dev\").",
       call. = FALSE)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

if (!file.exists("uSDT_pub.Rproj")) {
  stop("Run this script from the project root (the folder holding ",
       "uSDT_pub.Rproj).", call. = FALSE)
}

# theme_usdt() and save_figure()
source("R scripts/R functions/plot_functions.R")

dir.create("Figures/Manuscript", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 4: Settings
# ─────────────────────────────────────────────────────────────────────────────

seed <- 2025

# Direct sensitivity exceeds indirect sensitivity by .15.
# True individual differences limit power even with many trials.
population <- uSDT::usdt_population(
  dprime       = c(awareness = 0.30, influence = 0.15),  # mean true sensitivity
  sd_dprime    = 0.15,                                   # SD of true sensitivity
  cor_dprime   = 0,                                      # their correlation
  criterion    = 0,                                      # mean criterion
  sd_criterion = 0.30                                    # SD of criterion
)

# Trials per condition; usdt_power() expects totals for both conditions.
n_subjects           <- 40
trials_per_condition <- list(awareness = seq(100, 3, by = -1), influence = 150)

# Use "simulation" for full fitted studies; "approximate" is faster.
nsim    <- 1000
alpha   <- .05
method  <- "approximate"
n_cores <- 15

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 5: Power
# ─────────────────────────────────────────────────────────────────────────────

# Two-sided Wald test; nsim, seed and ncores only affect simulation.
power <- uSDT::usdt_power(population,
                          subjects   = n_subjects,
                          trials     = lapply(trials_per_condition,
                                              function(k) 2 * k),
                          hypothesis = "diff = 0",
                          alpha      = alpha,
                          method     = method,
                          nsim       = nsim,
                          seed       = seed,
                          ncores     = n_cores)

curve <- data.frame(trials_direct = power$trials_awareness / 2,
                    power         = power$power,
                    reliability   = power$reliability_awareness)
curve <- curve[order(curve$trials_direct), ]

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 6: Figure, power as a function of direct-task trials
# ─────────────────────────────────────────────────────────────────────────────

# Mark the first design reaching 80% power.
crossing <- curve[which(curve$power >= 0.80)[1], ]
crossing$label <- sprintf("%d trials", crossing$trials_direct)

fig <- ggplot(curve, aes(trials_direct, power)) +
  geom_hline(yintercept = 0.80, linetype = "13", colour = "grey50",
             linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(data = crossing, shape = 21, size = 1.7, fill = "white",
             stroke = 0.7) +
  geom_text(data = crossing,
            aes(x = trials_direct - 1, y = power + 0.07, label = label),
            hjust = 1, size = 2.5) +
  scale_x_continuous(breaks = unique(c(min(curve$trials_direct),
                                       curve$trials_direct[curve$trials_direct %% 10 == 0])),
                     expand = expansion(mult = c(0.02, 0.03))) +
  scale_y_continuous(labels = function(x) paste0(100 * x, "%"),
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Trials per condition in the awareness task",
       y = "Probability of a significant difference") +
  theme_usdt()

save_figure(fig, "Manuscript/Figure_5.png", width = 4.8, height = 3.1)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 7: Values reported in the manuscript
# ─────────────────────────────────────────────────────────────────────────────

cat("\n--- Power of d'(indirect) - d'(direct) ---\n")
cat(sprintf("  d'(direct) = %.2f, d'(indirect) = %.2f, N = %d, %d trials per condition in the indirect task\n",
            population$mu_A, population$mu_I, n_subjects,
            trials_per_condition$influence))
cat(sprintf("  %s power, alpha = %.2f; reliability of the indirect task = %.2f\n",
            method, alpha, power$reliability_influence[1]))
print(transform(curve, power = round(power, 3), reliability = round(reliability, 2)),
      row.names = FALSE)

cat(sprintf("\n  Fewest direct-task trials per condition with power >= .80: %d\n",
            crossing$trials_direct))

# ─────────────────────────────────────────────────────────────────────────────
