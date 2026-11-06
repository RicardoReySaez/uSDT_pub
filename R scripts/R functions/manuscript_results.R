# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : manuscript_results.R                                       ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.com                                 ║
# ║ Date          : 18-10-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Functions that read the results saved by the scripts in
# "R scripts/Empirical illustrations" and format them for the tables and the
# text of the manuscript, so that every reported value comes from the fitted
# models.

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Formatting
# ─────────────────────────────────────────────────────────────────────────────

# Number with a fixed number of decimals. lead = FALSE drops the leading zero,
# as APA style asks for values that cannot exceed 1 (correlations, p values).
fmt <- function(x, digits = 2, lead = TRUE) {
  x   <- round(x, digits)
  x[x == 0] <- 0
  out <- formatC(x, format = "f", digits = digits)
  if (!lead) out <- sub("^(-?)0\\.", "\\1.", out)
  out
}

# Interval in brackets
fmt_ci <- function(low, high, digits = 2, lead = TRUE) {
  sprintf("[%s, %s]", fmt(low, digits, lead), fmt(high, digits, lead))
}

# p value: "< .001" or three decimals without the leading zero. In the text,
# relation = TRUE adds "=" before the value, so that "p `r fmt_p(x, TRUE)`"
# reads "p = .016" or "p < .001".
fmt_p <- function(p, relation = FALSE) {
  ifelse(p < .001, "< .001",
         paste0(if (relation) "= " else "", fmt(p, 3, lead = FALSE)))
}

# Bayes factor for LaTeX math: "< .001" when tiny, two decimals when moderate,
# a thousands separator when large, and scientific notation when huge
fmt_bf <- function(x) {
  if (x < .001) return("< .001")
  if (x < 1000) return(fmt(x, 2))
  if (x < 1e5)  return(gsub(",", "{,}", format(round(x), big.mark = ",")))
  e <- floor(log10(x))
  sprintf("%s \\times 10^{%d}", fmt(x / 10^e, 1), e)
}

# Percentage of single trials classified correctly at a given d'
fmt_pc <- function(d) {
  sprintf("%.0f%%", 100 * pnorm(d / 2))
}

# Table cells: math minus signs, so that LaTeX does not print hyphens
tex <- function(x) {
  gsub("-", "$-$", x, fixed = TRUE)
}

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: Reading results
# ─────────────────────────────────────────────────────────────────────────────

# Results of one illustration. Paths are relative to the manuscript folder.
read_results <- function(name) {
  readRDS(file.path("..", "Results", "Empirical illustrations", paste0(name, ".rds")))
}

# Row of a uSDT test table (m$tests) for one hypothesis: H1 (difference), H2
# (correlation), or H3 (intercept)
test_row <- function(tests, hypothesis) {
  term <- c(H1 = "d'(influence) - d'(awareness)", H2 = "correlation",
            H3 = "intercept")[[hypothesis]]
  tests[tests$hypothesis == hypothesis & tests$term == term, ]
}

# Bayes factor row of usdt_bf() for one hypothesis, in the order H1, H2, H3
bf_row <- function(bf, hypothesis) {
  bf[match(hypothesis, c("H1", "H2", "H3")), ]
}

# Reliability of one task ("Awareness" or "Influence")
rel_row <- function(reliability, task) {
  reliability[reliability$task == task, ]
}

# ─────────────────────────────────────────────────────────────────────────────
