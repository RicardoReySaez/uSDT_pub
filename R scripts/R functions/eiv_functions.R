# ╔════════════════════════════════════════════════════════════════════════════╗
# ║                             SCRIPT OVERVIEW                                ║
# ╠════════════════════════════════════════════════════════════════════════════╣
# ║ Script Name   : eiv_functions.R                                            ║
# ║ Author        : Ricardo Rey-Sáez                                           ║
# ║ Role          : PhD Student in Psychology                                  ║
# ║ Institution   : Autonomous University of Madrid, Spain                     ║
# ║ Email         : ricardoreysaez95@gmail.es                                  ║
# ║ Date          : 02-11-2026                                                 ║
# ╚════════════════════════════════════════════════════════════════════════════╝

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 1: Description
# ─────────────────────────────────────────────────────────────────────────────

# Errors-in-variables regression functions (Klauer et al., 1998).
# Retrieved from: https://osf.io/36hk2/files/osfstorage (in the .zip folder!)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 2: Load Packages
# ─────────────────────────────────────────────────────────────────────────────

# minqa needed for estimating Hessian matrix
if (!requireNamespace("minqa", quietly = TRUE)) {
  install.packages("minqa")
}
library(minqa)

# ─────────────────────────────────────────────────────────────────────────────

# ─────────────────────────────────────────────────────────────────────────────
# SECTION 3: User-Defined Functions
# ─────────────────────────────────────────────────────────────────────────────

# Errors-in-variables log-likelihood function
eiv_loglik <- function(par, vis, cond) {
  mu <- par[1]; sigma <- par[2]; sx <- par[3]; sy <- par[4]
  beta <- par[5]; alpha <- par[6]
  nu    <- (sy * sigma)^2 + (beta * sx * sigma)^2 + (sx * sy)^2
  xi    <- (vis * sy / sx * sigma + beta * (cond - alpha) * sx / sy * sigma +
              mu * sx * sy / sigma) / sqrt(nu)
  gamma <- vis^2 / sx^2 + (cond - alpha)^2 / sy^2 + mu^2 / sigma^2
  f <- exp(-(gamma - xi^2) / 2) * pnorm(xi) / (2 * pi * sqrt(nu)) +
    exp(-(gamma - (mu / sigma)^2) / 2) * pnorm(-mu / sigma) / (2 * pi * sx * sy)
  ll <- -sum(log(f))
  if (!is.finite(ll)) 1e20 else ll
}

# Errors-in-variables estimate function
eiv_lm <- function(vis, cond, starts = 20) {
  fits <- lapply(seq_len(starts), function(i) {
    start <- runif(6, c(-1, 0.1, 0.1, 0.1, -0.2, -0.2), c(1, 2, 1, 1, 0.2, 0.2))
    minqa::bobyqa(start, eiv_loglik, vis = vis, cond = cond,
                  lower = c(-Inf, 1e-4, 1e-4, 1e-4, -Inf, -Inf))
  })
  best <- fits[[which.min(sapply(fits, `[[`, "fval"))]]
  se   <- sqrt(diag(solve(numDeriv::hessian(eiv_loglik, best$par,
                                            vis = vis, cond = cond))))
  data.frame(parameter = c("mu", "sigma", "sx", "sy", "slope", "intercept"),
             estimate  = best$par,
             se        = se)
}

# ─────────────────────────────────────────────────────────────────────────────
