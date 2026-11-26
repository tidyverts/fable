globalVariables(c("p", "q"))

#' @importFrom stats coef .lm.fit
specials_tbats <- new_specials(
  trend = function(method = c("N", "A", "Ad")) {
    if (!all(is.element(method, c("N", "A", "Ad")))) {
      cli::cli_abort("Invalid trend type")
    }
    list(method = method)
  },
  season = function(period = NULL, harmonics = NULL) {
    period <- get_frequencies(period, self$data, .auto = "all")
    period <- unname(period[period > 1])
    m <- length(period)
    if (is.null(harmonics)) {
      harmonics <- vector("list", m)
    } else if (!is.list(harmonics)) {
      harmonics <- as.list(rep_len(harmonics, m))
    } else if (length(harmonics) != m) {
      cli::cli_abort("`harmonics` must be a single value, or a list matching the number of seasonal `period`s")
    }
    list(period = period, harmonics = harmonics)
  },
  pdq = function(p = 0:5, q = 0:5) {
    list(p = p, q = q)
  },
  xreg = no_xreg,
  .required_specials = c("trend", "season", "pdq")
)

train_tbats <- function(.data, specials, ic = c("aicc", "aic", "bic"), ...) {
  ic <- match.arg(ic)
  ic_field <- c(aicc = "AICc", aic = "AIC", bic = "BIC")[[ic]]

  if (length(measured_vars(.data)) > 1) {
    cli::cli_abort("Only univariate responses are supported by TBATS.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]
  idx <- unclass(.data)[[index_var(.data)]]
  if (anyNA(y)) {
    cli::cli_abort("TBATS does not support missing values.")
  }

  trend_spec <- specials$trend[[1]]
  season_spec <- specials$season[[1]]
  pdq_spec <- specials$pdq[[1]]

  period <- season_spec$period
  # Harmonics are capped below the Nyquist frequency (k == period/2 exactly,
  # for even periods): at that harmonic the seasonal state's sine component
  # is identically zero (sin(pi) == 0), making the corresponding column of
  # the seed-state design matrix degenerate.
  max_k <- pmax(floor((period - 1) / 2), 1)
  harmonics <- map2(season_spec$harmonics, max_k, function(h, mx) {
    if (is.null(h)) {
      return(NULL)
    }
    kept <- h[h <= mx]
    if (length(kept) < length(h)) {
      cli::cli_warn("Reducing the number of Fourier harmonics to {mx} to avoid exceeding the Nyquist frequency.")
    }
    if (length(kept) == 0) mx else kept
  })

  # As with ARIMA()'s pdq(), avoid searching ARMA orders that can't be
  # estimated from the available data.
  max_pq <- max(floor(length(y) / 3), 0)
  pdq_spec$p <- pdq_spec$p[pdq_spec$p <= max_pq]
  pdq_spec$q <- pdq_spec$q[pdq_spec$q <= max_pq]

  trend_opts <- expand.grid(trend = trend_spec$method, stringsAsFactors = FALSE)
  trend_opts$use_trend <- trend_opts$trend != "N"
  trend_opts$use_damped <- trend_opts$trend == "Ad"

  candidates <- map2(trend_opts$use_trend, trend_opts$use_damped, function(use_trend, use_damped) {
    select_tbats_harmonics(y, use_trend, use_damped, period, harmonics)
  })
  ic_vals <- map_dbl(candidates, function(x) x$fit[[ic_field]])
  best <- candidates[[which.min(ic_vals)]]
  best_fit <- add_tbats_arma_errors(y, best$fit, p_range = pdq_spec$p, q_range = pdq_spec$q)

  par <- tbats_par_tibble(best_fit$par, best_fit$use_damped)
  state_names <- tbats_state_names(best_fit$use_trend, best$k, best_fit$p, best_fit$q)
  states_mat <- rbind(best_fit$x0, best_fit$states)
  colnames(states_mat) <- state_names

  structure(
    list(
      par = par,
      est = mutate(dplyr::ungroup(.data), .fitted = best_fit$fitted, .resid = best_fit$e),
      fit = tibble(
        sigma2 = best_fit$sigma2, log_lik = best_fit$log_lik,
        AIC = best_fit$AIC, AICc = best_fit$AICc, BIC = best_fit$BIC
      ),
      states = tsibble(
        !!!set_names(list(seq(idx[1] - index_step(idx, interval(.data)),
          by = index_step(idx, interval(.data)),
          length.out = NROW(states_mat)
        )), index_var(.data)),
        !!!set_names(split(states_mat, col(states_mat)), colnames(states_mat)),
        index = !!index(.data)
      ),
      spec = tibble(
        use_trend = best_fit$use_trend, damped = best_fit$use_damped,
        p = best_fit$p, q = best_fit$q,
        period = list(best$period), k = list(best$k)
      )
    ),
    class = "TBATS"
  )
}

# `par` (a nested list from tbats_unpack()/relist()) is flattened into a
# term/estimate tibble, giving each coefficient a readable name for
# tidy()/report() "for free" (rather than a separate labelling step).
tbats_par_tibble <- function(par, use_damped) {
  # alpha/beta/phi are always scalars, so their bare name is unambiguous;
  # gamma1/gamma2/ar/ma are always looked up by index (tbats_rebuild(),
  # tbats_get_par()), even when there's only one, so they're always indexed.
  scalar <- function(nm, x) if (is.null(x)) NULL else set_names(x, nm)
  indexed <- function(nm, x) if (is.null(x)) NULL else set_names(x, paste0(nm, seq_along(x)))
  vec <- c(
    scalar("alpha", par$alpha),
    scalar("beta", par$beta),
    if (use_damped) scalar("phi", par$phi),
    indexed("gamma1_", par$gamma1),
    indexed("gamma2_", par$gamma2),
    indexed("ar", par$ar),
    indexed("ma", par$ma)
  )
  tibble(term = names(vec), estimate = unname(vec))
}

tbats_state_names <- function(use_trend, k, p, q) {
  nm <- "l"
  if (use_trend) nm <- c(nm, "b")
  if (length(k)) {
    nm <- c(nm, unlist(lapply(seq_along(k), function(i) {
      c(paste0("season", i, "_", seq_len(k[i])), paste0("season", i, "_", seq_len(k[i]), "_star"))
    })))
  }
  if (p > 0) nm <- c(nm, paste0("ar", seq_len(p)))
  if (q > 0) nm <- c(nm, paste0("ma", seq_len(q)))
  nm
}

tbats_get_par <- function(object, nm) {
  object$par$estimate[object$par$term == nm]
}

# Rebuild the fitted state space matrices (and any parameters/dims implied by
# them) from a fitted TBATS object's `par`/`spec`, for use by forecast(),
# generate() and refit().
tbats_rebuild <- function(object) {
  spec <- object$spec
  k <- spec$k[[1]]
  period <- spec$period[[1]]
  vec <- function(nm, n) if (n == 0) NULL else vapply(seq_len(n), function(i) tbats_get_par(object, paste0(nm, i)), numeric(1))

  alpha <- tbats_get_par(object, "alpha")
  beta <- if (spec$use_trend) tbats_get_par(object, "beta") else NULL
  phi <- if (spec$use_trend) {
    if (spec$damped) tbats_get_par(object, "phi") else 1
  } else NULL
  gamma1 <- vec("gamma1_", length(k))
  gamma2 <- vec("gamma2_", length(k))
  ar <- vec("ar", spec$p)
  ma <- vec("ma", spec$q)

  ssm <- tbats_ssm(alpha, beta, phi, gamma1, gamma2, period, k, ar, ma)
  list(ssm = ssm, k = k, period = period, ar = ar, ma = ma)
}

#' TBATS model (Exponential smoothing state space model with Box-Cox
#' transformation, ARMA errors, Trend and Seasonal components)
#'
#' Fits a TBATS model applied to the (optionally transformed) response, as
#' described in De Livera, Hyndman & Snyder (2011). The fitted model is an
#' innovations state space model with a trigonometric representation of
#' seasonality, so a single `TBATS()` model covers what the original
#' methodology paper describes as both BATS and TBATS - a non-seasonal
#' specification, or one using the maximum number of harmonics for each
#' seasonal period (`season(harmonics = period/2)`), reproduce BATS'
#' seasonal-dummy behaviour exactly.
#'
#' Unlike the `forecast` package's `bats()`/`tbats()`, the Box-Cox
#' transformation is not estimated internally. Instead, as with other fable
#' models (e.g. `ETS(log(y) ~ ...)`), apply it via the model formula, for
#' example `TBATS(box_cox(y, lambda) ~ ...)`, with `lambda` chosen ahead of
#' time (e.g. using `feasts::guerrero()`). This lets fabletools' existing,
#' generic response-transformation handling (including bias-adjusted
#' back-transformed forecasts) apply, rather than duplicating it internally.
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used to select the trend
#' specification (from `trend()`) and the number of seasonal harmonics (from
#' `season()`, when not fixed).
#' @param ... Other arguments, currently ignored.
#'
#' @section Specials:
#'
#' \subsection{trend}{
#' The `trend` special is used to specify the form of the trend term.
#' \preformatted{
#' trend(method = c("N", "A", "Ad"))
#' }
#' \tabular{ll}{
#'   `method` \tab The form of the trend: none ("N"), additive ("A") or
#'   damped additive ("Ad"). All specified methods are fitted, and the one
#'   minimising `ic` is kept.
#' }
#' }
#'
#' \subsection{season}{
#' The `season` special specifies the seasonal periods and the number of
#' trigonometric (Fourier) harmonics used to represent each of them.
#' \preformatted{
#' season(period = NULL, harmonics = NULL)
#' }
#' \tabular{ll}{
#'   `period`    \tab The periodic nature of the seasonality(s). This can be
#'   a vector, allowing multiple seasonal periods (as in `fourier()`). If
#'   omitted, all seasonal periods implied by the data's interval are used.\cr
#'   `harmonics` \tab The number of Fourier terms used for each seasonal
#'   `period`. If omitted (the default), the number of harmonics per period
#'   is chosen automatically by minimising `ic`.
#' }
#' }
#'
#' \subsection{pdq}{
#' The `pdq` special specifies the search range for the order of the ARMA(p,q)
#' process used to model the errors.
#' \preformatted{
#' pdq(p = 0:5, q = 0:5)
#' }
#' \tabular{ll}{
#'   `p` \tab The candidate autoregressive order(s) to search over.\cr
#'   `q` \tab The candidate moving average order(s) to search over.
#' }
#' Set `pdq(p = 0, q = 0)` to disable ARMA errors.
#' }
#'
#' @return A model specification.
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value))
#' @seealso
#' [Forecasting: Principles and Practices, Complex seasonality (chapter 12)](https://otexts.com/fpp3/complexseasonality.html)
#'
#' @references De Livera, A.M., Hyndman, R.J., & Snyder, R. D. (2011),
#' Forecasting time series with complex seasonal patterns using exponential
#' smoothing, \emph{Journal of the American Statistical Association},
#' \bold{106}(496), 1513-1527.
#'
#' @export
TBATS <- function(formula, ic = c("aicc", "aic", "bic"), ...) {
  ic <- match.arg(ic)
  tbats_model <- new_model_class("TBATS",
    train = train_tbats, specials = specials_tbats,
    check = all_tsbl_checks
  )
  new_model_definition(tbats_model, !!enquo(formula), ic = ic, ...)
}

#' @inherit forecast.ARIMA
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value)) %>%
#'   forecast()
#' @export
forecast.TBATS <- function(object, new_data, specials = NULL, ...) {
  h <- NROW(new_data)
  rb <- tbats_rebuild(object)
  F <- rb$ssm$F
  g <- rb$ssm$g
  w <- rb$ssm$w

  laststate <- as.numeric(object$states[NROW(object$states), tbats_state_names(object$spec$use_trend, rb$k, object$spec$p, object$spec$q)])

  x <- laststate
  yhat <- numeric(h)
  for (t in seq_len(h)) {
    yhat[t] <- sum(w * x)
    x <- as.numeric(F %*% x)
  }

  vm <- numeric(h)
  vm[1] <- 1
  if (h > 1) {
    Fg <- g # F^(j-1) g, starting at j = 1
    for (j in seq_len(h - 1)) {
      cj <- sum(w * Fg)
      vm[j + 1] <- vm[j] + cj^2
      Fg <- as.numeric(F %*% Fg)
    }
  }

  distributional::dist_normal(yhat, sqrt(object$fit$sigma2 * vm))
}

#' @inherit generate.ETS
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value)) %>%
#'   generate(times = 100)
#' @export
generate.TBATS <- function(x, new_data, specials, bootstrap = FALSE, ...) {
  if (!is_regular(new_data)) {
    cli::cli_abort("Simulation new_data must be regularly spaced")
  }
  start_idx <- min(new_data[[index_var(new_data)]])
  start_pos <- match(start_idx - index_step(start_idx, interval(new_data)), x$states[[index_var(x$states)]])
  if (is.na(start_pos)) {
    cli::cli_abort("The first observation index of simulation data must be within the model's training set.")
  }

  rb <- tbats_rebuild(x)
  state_names <- tbats_state_names(x$spec$use_trend, rb$k, x$spec$p, x$spec$q)
  initstate <- as.numeric(x$states[start_pos, state_names])

  if (!(".innov" %in% names(new_data))) {
    if (bootstrap) {
      new_data$.innov <- sample(stats::na.omit(residuals(x) - mean(residuals(x), na.rm = TRUE)),
        NROW(new_data),
        replace = TRUE
      )
    } else {
      new_data$.innov <- stats::rnorm(NROW(new_data), sd = sqrt(x$fit$sigma2))
    }
  }

  F <- rb$ssm$F
  g <- rb$ssm$g
  w <- rb$ssm$w
  sim_path <- function(innov) {
    xp <- initstate
    out <- numeric(length(innov))
    for (t in seq_along(out)) {
      out[t] <- sum(w * xp) + innov[t]
      xp <- as.numeric(F %*% xp) + g * innov[t]
    }
    out
  }

  transmute(group_by_key(new_data), ".sim" := sim_path(!!sym(".innov")))
}

#' Refit a TBATS model
#'
#' Applies a fitted TBATS model to a new dataset.
#'
#' @inheritParams refit.ARIMA
#' @param reinitialise If TRUE, all parameters (including the seasonal
#' harmonics and ARMA error order) will be re-estimated to suit the new data.
#' If FALSE, the fitted structure and parameters are kept, and only applied
#' (filtered) to the new data.
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(TBATS(value))
#'
#' fit %>%
#'   refit(lung_deaths_female) %>%
#'   report()
#' @export
refit.TBATS <- function(object, new_data, specials = NULL, reinitialise = TRUE, ...) {
  y <- unclass(new_data)[[measured_vars(new_data)]]
  idx <- unclass(new_data)[[index_var(new_data)]]

  if (reinitialise) {
    spec <- object$spec
    sel <- select_tbats_harmonics(y, spec$use_trend, spec$damped, spec$period[[1]],
      harmonics = as.list(spec$k[[1]])
    )
    best_fit <- if (spec$p > 0 || spec$q > 0) {
      estimate_tbats(y, spec$use_trend, spec$damped, sel$period, sel$k,
        ar = numeric(spec$p), ma = numeric(spec$q))
    } else sel$fit

    par <- tbats_par_tibble(best_fit$par, best_fit$use_damped)
    state_names <- tbats_state_names(best_fit$use_trend, sel$k, best_fit$p, best_fit$q)
    states_mat <- rbind(best_fit$x0, best_fit$states)
    colnames(states_mat) <- state_names

    object$par <- par
    object$fit <- tibble(
      sigma2 = best_fit$sigma2, log_lik = best_fit$log_lik,
      AIC = best_fit$AIC, AICc = best_fit$AICc, BIC = best_fit$BIC
    )
    object$spec <- tibble(
      use_trend = best_fit$use_trend, damped = best_fit$use_damped,
      p = best_fit$p, q = best_fit$q, period = list(sel$period), k = list(sel$k)
    )
    object$est <- mutate(dplyr::ungroup(new_data), .fitted = best_fit$fitted, .resid = best_fit$e)
    object$states <- tsibble(
      !!!set_names(list(seq(idx[1] - index_step(idx, interval(new_data)),
        by = index_step(idx, interval(new_data)),
        length.out = NROW(states_mat)
      )), index_var(new_data)),
      !!!set_names(split(states_mat, col(states_mat)), colnames(states_mat)),
      index = !!index(new_data)
    )
    return(object)
  }

  rb <- tbats_rebuild(object)
  state_names <- tbats_state_names(object$spec$use_trend, rb$k, object$spec$p, object$spec$q)
  x0 <- as.numeric(object$states[1, state_names])
  fit <- tbats_filter(y, rb$ssm$F, rb$ssm$g, rb$ssm$w, x0)
  n <- length(y)
  n_par <- NROW(object$par) + length(rb$ssm$w)

  states_mat <- rbind(x0, fit$states)
  colnames(states_mat) <- state_names

  object$est <- mutate(dplyr::ungroup(new_data), .fitted = fit$fitted, .resid = fit$e)
  object$states <- tsibble(
    !!!set_names(list(seq(idx[1] - index_step(idx, interval(new_data)),
      by = index_step(idx, interval(new_data)),
      length.out = NROW(states_mat)
    )), index_var(new_data)),
    !!!set_names(split(states_mat, col(states_mat)), colnames(states_mat)),
    index = !!index(new_data)
  )
  lik <- n * log(sum(fit$e^2))
  object$fit <- tibble(
    sigma2 = sum(fit$e^2) / n, log_lik = -0.5 * (lik + n * (log(2 * pi / n) + 1)),
    AIC = lik + 2 * n_par, BIC = lik + log(n) * n_par,
    AICc = lik + 2 * n_par + 2 * n_par * (n_par + 1) / (n - n_par - 1)
  )
  object
}

#' Extend a fitted TBATS model with new data
#'
#' Applies a fitted TBATS model to a new, extended dataset. `new_data` must
#' start immediately after the end of `object`'s training data. The model's
#' structure and coefficients are held fixed - only the fitted values, states
#' and fit statistics are updated to account for the new observations.
#'
#' @inheritParams stream.ARIMA
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   filter(index < yearmonth("1979 Jan")) %>%
#'   model(TBATS(value))
#'
#' fit %>%
#'   stream(lung_deaths_male %>% filter(index >= yearmonth("1979 Jan"))) %>%
#'   report()
#' @export
stream.TBATS <- function(object, new_data, specials = NULL, ...) {
  idx_var <- index_var(object$est)
  n_old <- NROW(object$est)
  stream_start <- object$est[[idx_var]][n_old] + index_step(object$est[[idx_var]], interval(object$est))
  if (unclass(new_data)[[index_var(new_data)]][1] != stream_start) {
    cli::cli_abort("Streaming to a TBATS model must start one step beyond the end of the trained data.")
  }

  y_new <- unclass(new_data)[[measured_vars(new_data)]]
  if (anyNA(y_new)) {
    cli::cli_abort("TBATS does not support missing values.")
  }

  rb <- tbats_rebuild(object)
  state_names <- tbats_state_names(object$spec$use_trend, rb$k, object$spec$p, object$spec$q)
  x0 <- as.numeric(object$states[NROW(object$states), state_names])

  fit <- tbats_filter(y_new, rb$ssm$F, rb$ssm$g, rb$ssm$w, x0)

  n <- n_old + length(y_new)
  n_par <- NROW(object$par) + length(rb$ssm$w)
  sse <- object$fit$sigma2 * n_old + sum(fit$e^2)
  lik <- n * log(sse)

  new_states <- fit$states
  colnames(new_states) <- state_names

  object$est <- dplyr::bind_rows(
    object$est,
    mutate(new_data, .fitted = fit$fitted, .resid = fit$e)
  )
  object$states <- dplyr::bind_rows(
    object$states,
    as_tibble(new_states) %>% mutate(!!idx_var := new_data[[index_var(new_data)]])
  )
  object$fit <- tibble(
    sigma2 = sse / n, log_lik = -0.5 * (lik + n * (log(2 * pi / n) + 1)),
    AIC = lik + 2 * n_par, BIC = lik + log(n) * n_par,
    AICc = lik + 2 * n_par + 2 * n_par * (n_par + 1) / (n - n_par - 1)
  )
  object
}

#' @inherit fitted.ARIMA
#' @export
fitted.TBATS <- function(object, ...) {
  object$est[[".fitted"]]
}

#' @inherit residuals.ARIMA
#' @export
residuals.TBATS <- function(object, ...) {
  object$est[[".resid"]]
}

#' Glance a TBATS model
#'
#' Construct a single row summary of the TBATS model.
#'
#' Contains the variance of residuals (`sigma2`), the log-likelihood
#' (`log_lik`), and information criteria (`AIC`, `AICc`, `BIC`).
#'
#' @inheritParams generics::glance
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value)) %>%
#'   glance()
#' @export
glance.TBATS <- function(x, ...) {
  x$fit
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value)) %>%
#'   tidy()
#' @export
tidy.TBATS <- function(x, ...) {
  x$par
}

#' @export
model_sum.TBATS <- function(x) {
  spec <- x$spec
  period <- spec$period[[1]]
  k <- spec$k[[1]]
  seas <- if (length(period)) {
    paste0("{", paste(sprintf("<%s,%s>", round(period, 2), k), collapse = ", "), "}")
  } else "{-}"
  phi <- if (!spec$use_trend) {
    "-"
  } else if (spec$damped) {
    format(round(tbats_get_par(x, "phi"), 3))
  } else {
    "1"
  }
  sprintf("TBATS({%s,%s}, %s, %s)", spec$p, spec$q, phi, seas)
}

#' @export
report.TBATS <- function(object, ...) {
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", format(tbats_get_par(object, "alpha")), "\n"))
  if (object$spec$use_trend) {
    cat(paste("    beta  =", format(tbats_get_par(object, "beta")), "\n"))
  }
  k <- object$spec$k[[1]]
  if (length(k)) {
    for (i in seq_along(k)) {
      cat(paste0("    gamma1_", i, " = ", format(tbats_get_par(object, paste0("gamma1_", i))), "\n"))
      cat(paste0("    gamma2_", i, " = ", format(tbats_get_par(object, paste0("gamma2_", i))), "\n"))
    }
  }
  if (object$spec$damped) {
    cat(paste("    phi   =", format(tbats_get_par(object, "phi")), "\n"))
  }
  if (object$spec$p > 0 || object$spec$q > 0) {
    cat(paste0("\n  ARMA(", object$spec$p, ",", object$spec$q, ") errors\n"))
  }

  cat("\n  sigma^2:  ")
  cat(round(object$fit$sigma2, 4))
  cat("\n\n")
  print(c(AIC = object$fit$AIC, AICc = object$fit$AICc, BIC = object$fit$BIC))
}

#' Extract estimated states from a TBATS model.
#'
#' @param object An estimated model.
#' @param ... Unused.
#'
#' @return A [fabletools::dable()] containing estimated states.
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(tbats = TBATS(value)) %>%
#'   components()
#' @export
components.TBATS <- function(object, ...) {
  spec <- object$spec
  period <- spec$period[[1]]
  k <- spec$k[[1]]
  response <- measured_vars(object$est)[[1]]

  out <- object$states[, c(index_var(object$states), "l", if (spec$use_trend) "b")]
  colnames(out)[colnames(out) == "l"] <- "level"
  if (spec$use_trend) colnames(out)[colnames(out) == "b"] <- "slope"

  seasonalities <- list()
  if (length(period)) {
    for (i in seq_along(period)) {
      cols <- paste0("season", i, "_", seq_len(k[i]))
      out[[paste0("season", i)]] <- rowSums(as.matrix(object$states[cols]))
    }
    seasonalities <- set_names(
      lapply(period, function(m) list(period = m, base = NA_real_)),
      paste0("season", seq_along(period))
    )
  }

  est_vars <- transmute(object$est, !!sym(response), remainder = !!sym(".resid"))
  out <- left_join(out, est_vars, by = index_var(object$states))

  # Each state's row is the *post-update* state (following ETS's convention),
  # so reconstructing the response at time t needs the previous row's states:
  # y[t] = level[t-1] (+ trend[t-1]) + sum(season_i[t-1]) + remainder[t].
  eqn <- expr(lag(!!sym("level"), 1))
  if (spec$use_trend) {
    phi <- if (spec$damped) tbats_get_par(object, "phi") else 1
    eqn <- expr(!!eqn + !!phi * lag(!!sym("slope"), 1))
  }
  for (i in seq_along(period)) {
    eqn <- expr(!!eqn + lag(!!sym(paste0("season", i)), 1))
  }
  eqn <- expr(!!eqn + !!sym("remainder"))

  fabletools::as_dable(out,
    resp = !!sym(response), method = model_sum(object),
    seasons = seasonalities, aliases = list2(!!response := eqn)
  )
}
