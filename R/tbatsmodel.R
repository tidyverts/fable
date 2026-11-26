# TBATS state space model machinery.
#
# TBATS represents the (possibly Box-Cox transformed) observation as a linear
# innovations state space model:
#   y[t] = w'x[t-1] + e[t]
#   x[t] = F x[t-1] + g e[t]
# with state x = (level, slope?, <trig seasonal terms>, <AR errors>, <MA errors>).
# BATS is the special case where each seasonal period's trigonometric terms
# use the maximum number of harmonics (k = floor(period/2)), so a single
# implementation of this state space covers both models - there is no need
# for a separate BATS code path. Any Box-Cox transformation of the response
# is handled by fabletools via the model formula (see R/tbats.R); everything
# here operates on `y` exactly as given.

# Build the F, g, w matrices for a TBATS state space model from its
# structural parameters. Rather than forecast's makeTBATSFMatrix(), which
# assembles these via a long combinatorial chain of cbind()/rbind() calls (one
# branch per combination of trend/seasonal/AR/MA being present or absent),
# this fixes the index ranges for each block up front and assigns directly
# into pre-allocated zero matrices/vectors.
tbats_ssm <- function(alpha, beta = NULL, phi = NULL, gamma1 = NULL, gamma2 = NULL,
                      period = NULL, k = NULL, ar = NULL, ma = NULL) {
  adj_beta <- !is.null(beta)
  tau <- if (!is.null(k)) 2L * sum(k) else 0L
  p <- length(ar)
  q <- length(ma)
  m <- 1L + adj_beta + tau + p + q

  w <- numeric(m)
  g <- numeric(m)
  F <- matrix(0, m, m)

  w[1] <- 1
  g[1] <- alpha
  F[1, 1] <- 1
  if (adj_beta) {
    w[2] <- phi
    g[2] <- beta
    F[1, 2] <- phi
    F[2, 2] <- phi
  }

  seas_idx <- integer()
  gamma_bold <- numeric(tau)
  if (tau > 0) {
    seas_idx <- (1L + adj_beta) + seq_len(tau)
    A <- matrix(0, tau, tau)
    pos <- 0L
    for (i in seq_along(k)) {
      idx <- pos + seq_len(k[i])
      w[seas_idx[idx]] <- 1
      gamma_bold[idx] <- gamma1[i]
      gamma_bold[idx + k[i]] <- gamma2[i]

      theta <- 2 * pi * seq_len(k[i]) / period[i]
      C <- diag(cos(theta), k[i])
      S <- diag(sin(theta), k[i])
      A[idx, idx] <- C
      A[idx, idx + k[i]] <- S
      A[idx + k[i], idx] <- -S
      A[idx + k[i], idx + k[i]] <- C

      pos <- pos + 2L * k[i]
    }
    g[seas_idx] <- gamma_bold
    F[seas_idx, seas_idx] <- A
  }

  if (p > 0) {
    ar_idx <- 1L + adj_beta + tau + seq_len(p)
    w[ar_idx] <- ar
    g[ar_idx[1]] <- 1
    F[1, ar_idx] <- alpha * ar
    if (adj_beta) F[2, ar_idx] <- beta * ar
    if (tau > 0) F[seas_idx, ar_idx] <- outer(gamma_bold, ar)
    F[ar_idx[1], ar_idx] <- ar
    if (p > 1) F[ar_idx[-1], ar_idx[-p]] <- diag(p - 1)
  }

  if (q > 0) {
    ma_idx <- 1L + adj_beta + tau + p + seq_len(q)
    w[ma_idx] <- ma
    g[ma_idx[1]] <- 1
    F[1, ma_idx] <- alpha * ma
    if (adj_beta) F[2, ma_idx] <- beta * ma
    if (tau > 0) F[seas_idx, ma_idx] <- outer(gamma_bold, ma)
    if (p > 0) F[ar_idx[1], ma_idx] <- ma
    if (q > 1) F[ma_idx[-1], ma_idx[-q]] <- diag(q - 1)
  }

  list(F = F, g = g, w = w)
}

# Roots "at infinity" (all-zero coefficients, e.g. the optimiser's starting
# point) leave `polyroot()` nothing to return - an empty root set is
# trivially outside the unit circle, so it's treated as admissible rather
# than triggering min()'s "no non-missing arguments" warning.
root_mod_ok <- function(coef, sign = 1, tol = 1e-2) {
  if (!length(coef)) {
    return(TRUE)
  }
  roots <- polyroot(c(1, sign * coef))
  !length(roots) || min(Mod(roots)) > 1 + tol
}

# The spectral radius of D = F - g w' must stay below 1 for the model to be
# stable, and AR/MA polynomials must be stationary/invertible. This single
# eigenvalue check subsumes the separate alpha/beta/gamma bound checks that
# would otherwise be needed, mirroring the "admissible" parameter region
# used for ETS models.
tbats_admissible <- function(F, g, w, ar = NULL, ma = NULL, tol = 1e-2) {
  if (!root_mod_ok(ar, -1, tol) || !root_mod_ok(ma, 1, tol)) {
    return(FALSE)
  }
  D <- F - outer(g, w)
  max(Mod(eigen(D, symmetric = FALSE, only.values = TRUE)$values)) < 1 + tol
}

# Parameters are packed into/out of a plain numeric vector using a named-list
# "template" (via utils::relist()), rather than hand-tracked integer offsets
# into a position-dependent vector plus a side-channel "control" list (as
# forecast's parameterise()/unParameteriseTBATS() do). The damping parameter
# phi is optimised on a logit scale so that it always maps back into its
# admissible range (0.8, 1) without needing a separate bounds check or
# rejection (returning Inf) during optimisation.
tbats_par_template <- function(use_trend, use_damped, k, p, q) {
  template <- list(alpha = 0.09)
  if (use_trend) {
    template$beta <- 0.05
    if (use_damped) template$phi <- stats::qlogis((0.999 - 0.8) / 0.2)
  }
  if (length(k)) {
    template$gamma1 <- rep(1e-3, length(k))
    template$gamma2 <- rep(1e-3, length(k))
  }
  if (p > 0) template$ar <- rep(0, p)
  if (q > 0) template$ma <- rep(0, q)
  template
}

tbats_unpack <- function(par, template) {
  par <- utils::relist(par, template)
  if (!is.null(par$phi)) {
    par$phi <- 0.8 + 0.2 * stats::plogis(par$phi)
  } else if (!is.null(par$beta)) {
    par$phi <- 1 # trend used but undamped
  }
  par
}

# Conditional on a set of structural (hyper)parameters, the seed state x0
# that minimises the sum of squared innovations is available in closed form:
# starting from x0 = 0, the effect of a nonzero x0 on e[t] is -w'D^(t-1) x0
# (D = F - gw'), so regressing the x0 = 0 residuals onto w'D^(t-1) recovers
# the optimal x0 by OLS. This keeps x0 out of the numerical optimisation
# entirely (computed once, held fixed while alpha/beta/gamma/ar/ma are
# optimised), similar in spirit to how ARIMA concentrates regression
# coefficients out of the likelihood.
tbats_seed_states <- function(y, F, g, w) {
  m <- length(w)
  n <- length(y)
  e0 <- tbats_filter(y, F, g, w, numeric(m))$e

  D <- F - outer(g, w)
  design <- matrix(0, n, m)
  row <- w
  for (t in seq_len(n)) {
    design[t, ] <- row
    row <- as.numeric(row %*% D)
  }
  as.numeric(.lm.fit(design, e0)$coefficients)
}

tbats_loglik <- function(par, template, y, period, k, x0) {
  par <- tbats_unpack(par, template)
  ssm <- tbats_ssm(
    alpha = par$alpha, beta = par$beta, phi = par$phi,
    gamma1 = par$gamma1, gamma2 = par$gamma2,
    period = period, k = k, ar = par$ar, ma = par$ma
  )
  if (!tbats_admissible(ssm$F, ssm$g, ssm$w, par$ar, par$ma)) {
    return(Inf)
  }
  e <- tbats_filter(y, ssm$F, ssm$g, ssm$w, x0)$e
  ll <- length(y) * log(sum(e^2))
  if (!is.finite(ll)) Inf else ll
}

# Fit a single, fully specified TBATS structure to `y`. `y` should already
# have any response transformation applied by the caller. Returns the fitted
# structural parameters, seed states, fitted values/residuals and IC (used to
# compare candidate structures).
estimate_tbats <- function(y, use_trend, use_damped, period = NULL, k = NULL,
                           ar = NULL, ma = NULL) {
  p <- length(ar)
  q <- length(ma)
  template <- tbats_par_template(use_trend, use_damped, k, p, q)
  par <- unlist(template, use.names = FALSE)

  beta0 <- if (use_trend) 0 else NULL
  phi0 <- if (use_trend && use_damped) 0.999 else if (use_trend) 1 else NULL
  gamma10 <- if (length(k)) rep(1e-3, length(k)) else NULL
  gamma20 <- if (length(k)) rep(1e-3, length(k)) else NULL
  ssm0 <- tbats_ssm(0.09, beta0, phi0, gamma10, gamma20, period, k, numeric(p), numeric(q))
  x0 <- tbats_seed_states(y, ssm0$F, ssm0$g, ssm0$w)

  opt <- stats::optim(
    par, tbats_loglik, method = "Nelder-Mead",
    template = template, y = y, period = period, k = k, x0 = x0,
    control = list(maxit = 500 * length(par))
  )
  loglik <- opt$value
  par <- tbats_unpack(opt$par, template)

  ssm <- tbats_ssm(par$alpha, par$beta, par$phi, par$gamma1, par$gamma2, period, k, par$ar, par$ma)
  fit <- tbats_filter(y, ssm$F, ssm$g, ssm$w, x0)

  n <- length(y)
  n_par <- length(x0) + length(par)
  aic <- loglik + 2 * n_par
  aicc <- aic + 2 * n_par * (n_par + 1) / (n - n_par - 1)
  bic <- loglik + log(n) * n_par

  list(
    par = par, x0 = x0, period = period, k = k, p = p, q = q,
    use_trend = use_trend, use_damped = use_damped,
    ar = par$ar, ma = par$ma,
    fitted = fit$fitted, e = fit$e, states = fit$states,
    sigma2 = sum(fit$e^2) / n,
    log_lik = -0.5 * (loglik + n * (log(2 * pi / n) + 1)),
    AIC = aic, AICc = aicc, BIC = bic, n_par = n_par
  )
}

# Choose the number of Fourier harmonics for each seasonal period by a
# coordinate-wise AICc scan (one period at a time, in increasing order of
# period), capping each period's harmonics so that two periods' terms can't
# represent the same seasonal frequency (which would make the model
# unidentifiable) - replacing forecast's three-branch up/level/down stepping
# search with a single monotone scan per period, stopping as soon as AICc
# stops improving.
select_tbats_harmonics <- function(y, use_trend, use_damped, period, harmonics,
                                   ar = NULL, ma = NULL, max_k = 6) {
  ord <- order(period)
  period <- period[ord]
  harmonics <- harmonics[ord]

  k <- vapply(harmonics, function(h) if (length(h) == 1) as.integer(h) else 1L, integer(1))
  auto <- vapply(harmonics, length, integer(1)) != 1

  best_fit <- NULL
  for (i in which(auto)) {
    max_ki <- min(max_k, floor((period[i] - 1) / 2))
    cand <- harmonics[[i]] %||% seq_len(max(max_ki, 1))
    cand <- cand[cand <= max_ki]
    earlier <- period[seq_len(i - 1)]
    cand <- Filter(function(kk) {
      period[i] %% kk != 0 || !any(earlier %% (period[i] / kk) == 0)
    }, cand)

    best_k <- k[i]
    for (kk in cand) {
      k[i] <- kk
      fit <- try(estimate_tbats(y, use_trend, use_damped, period, k, ar, ma), silent = TRUE)
      if (inherits(fit, "try-error")) next
      if (is.null(best_fit) || fit$AICc < best_fit$AICc) {
        best_fit <- fit
        best_k <- kk
      } else if (kk > best_k) {
        break
      }
    }
    k[i] <- best_k
  }

  if (is.null(best_fit)) {
    best_fit <- estimate_tbats(y, use_trend, use_damped, period, k, ar, ma)
  }
  list(k = k, period = period, fit = best_fit)
}

# Refine a fitted structure by testing for ARMA(p,q) errors in its residuals,
# mirroring forecast::tbats()'s use of an automatic ARMA order search as a
# final refinement step (rather than as part of the core structural search
# over trend/damping/harmonics). `stats::arima()` (already used by fable's
# `ARIMA()`) replaces forecast's dependency on `auto.arima()`.
add_tbats_arma_errors <- function(y, fit, p_range = 0:5, q_range = 0:5, order_constraint = 5) {
  model_opts <- expand.grid(p = p_range, q = q_range)
  model_opts <- model_opts[model_opts$p + model_opts$q > 0 & model_opts$p + model_opts$q <= order_constraint, ]

  best_ic <- Inf
  best_coef <- NULL
  best_pq <- c(0L, 0L)
  for (i in seq_len(NROW(model_opts))) {
    p <- model_opts$p[i]
    q <- model_opts$q[i]
    arma <- try(
      stats::arima(fit$e, order = c(p, 0, q), include.mean = FALSE, method = "ML"),
      silent = TRUE
    )
    if (inherits(arma, "try-error") || !is.finite(arma$aic)) next
    if (arma$aic < best_ic) {
      best_ic <- arma$aic
      best_pq <- c(p, q)
      best_coef <- stats::coef(arma)
    }
  }
  if (!is.finite(best_ic)) {
    return(fit)
  }

  ar0 <- if (best_pq[1] > 0) unname(best_coef[seq_len(best_pq[1])]) else NULL
  ma0 <- if (best_pq[2] > 0) unname(best_coef[best_pq[1] + seq_len(best_pq[2])]) else NULL

  candidate <- try(
    estimate_tbats(y, fit$use_trend, fit$use_damped, fit$period, fit$k, ar0, ma0),
    silent = TRUE
  )
  if (inherits(candidate, "try-error") || candidate$AICc >= fit$AICc) fit else candidate
}
