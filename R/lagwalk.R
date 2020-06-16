train_lagwalk <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by lagwalks.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]
  n <- length(y)

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }

  drift <- specials$drift[[1]] %||% FALSE
  lag <- specials$lag[[1]]

  y_na <- which(is.na(y))
  y_na <- y_na[y_na > lag]
  fits <- stats::lag(y, -lag)
  for (i in y_na) {
    if (is.na(fits)[i]) {
      fits[i] <- fits[i - lag]
    }
  }

  fitted <- c(rep(NA, min(lag, n)), utils::head(fits, -lag))
  if (drift) {
    fit <- summary(stats::lm(y - fitted ~ 1, na.action = stats::na.exclude))
    b <- fit$coefficients[1, 1]
    b.se <- fit$coefficients[1, 2]
    sigma <- fit$sigma
    fitted <- fitted + b
  }
  else {
    b <- b.se <- dbl()
    sigma <- stats::sd(y - fitted, na.rm = TRUE)
  }
  res <- y - fitted

  structure(
    list(
      b = b,
      b.se = b.se,
      lag = lag,
      sigma2 = sigma^2,
      .fitted = fitted,
      .resid = res,
      time = list(start = unclass(.data)[[index_var(.data)]][[1]], interval = interval(.data)),
      future = y[c(rep(NA, max(0, lag - n)), seq_len(min(n, lag)) + n - min(n, lag))]
    ),
    class = "RW"
  )
}

#' Random walk models
#'
#' \code{RW()} returns a random walk model, which is equivalent to an ARIMA(0,1,0)
#' model with an optional drift coefficient included using \code{drift()}. \code{naive()} is simply a wrapper
#' to \code{rwf()} for simplicity. \code{snaive()} returns forecasts and
#' prediction intervals from an ARIMA(0,0,0)(0,1,0)m model where m is the
#' seasonal period.
#'
#' The random walk with drift model is \deqn{Y_t=c + Y_{t-1} + Z_t}{Y[t]=c +
#' Y[t-1] + Z[t]} where \eqn{Z_t}{Z[t]} is a normal iid error. Forecasts are
#' given by \deqn{Y_n(h)=ch+Y_n}{Y[n+h]=ch+Y[n]}. If there is no drift (as in
#' \code{naive}), the drift parameter c=0. Forecast standard errors allow for
#' uncertainty in estimating the drift parameter (unlike the corresponding
#' forecasts obtained by fitting an ARIMA model directly).
#'
#' The seasonal naive model is \deqn{Y_t= Y_{t-m} + Z_t}{Y[t]=Y[t-m] + Z[t]}
#' where \eqn{Z_t}{Z[t]} is a normal iid error.
#'
#' @aliases report.RW
#'
#' @param formula Model specification (see "Specials" section).
#' @param ... Not used.
#'
#' @section Specials:
#'
#' \subsection{lag}{
#' The `lag` special is used to specify the lag order for the random walk process.
#' If left out, this special will automatically be included.
#'
#' \preformatted{
#' lag(lag = NULL)
#' }
#'
#' \tabular{ll}{
#'   `lag`        \tab The lag order for the random walk process. If `lag = m`, forecasts will return the observation from `m` time periods ago. This can also be provided as text indicating the duration of the lag window (for example, annual seasonal lags would be "1 year").
#' }
#' }
#'
#' \subsection{drift}{
#' The `drift` special can be used to include a drift/trend component into the model. By default, drift is not included unless `drift()` is included in the formula.
#'
#' \preformatted{
#' drift(drift = TRUE)
#' }
#'
#' \tabular{ll}{
#'   `drift`      \tab If `drift = TRUE`, a drift term will be included in the model.
#' }
#' }
#'
#' @return A model specification.
#'
#' @seealso
#' [Forecasting: Principles and Practices, Some simple forecasting methods (section 3.2)](https://otexts.com/fpp3/simple-methods.html)
#'
#' @examples
#' library(tsibbledata)
#' aus_production %>%
#'   model(rw = RW(Beer ~ drift()))
#' @export
RW <- function(formula, ...) {
  rw_model <- new_model_class("RW",
    train = train_lagwalk,
    specials = new_specials(
      lag = function(lag = NULL) {
        if (is.null(lag)) {
          lag <- 1
        }
        if (!rlang::is_integerish(lag)) {
          warn("Non-integer lag orders for random walk models are not supported. Rounding to the nearest integer.")
          lag <- round(lag)
        }
        get_frequencies(lag, self$data, .auto = "smallest")
      },
      drift = function(drift = TRUE) {
        drift
      },
      xreg = no_xreg,
      .required_specials = c("lag")
    ),
    check = all_tsbl_checks
  )
  new_model_definition(rw_model, !!enquo(formula), ...)
}

#' @rdname RW
#'
#' @examples
#'
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value))
#' @export
NAIVE <- RW

#' @rdname RW
#'
#' @examples
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year")))
#' @export
SNAIVE <- function(formula, ...) {
  snaive_model <- new_model_class("RW",
    train = train_lagwalk,
    specials = new_specials(
      lag = function(lag = NULL) {
        lag <- get_frequencies(lag, self$data, .auto = "smallest")
        if (lag == 1) {
          abort("Non-seasonal model specification provided, use RW() or provide a different lag specification.")
        }
        if (!rlang::is_integerish(lag)) {
          warn("Non-integer lag orders for random walk models are not supported. Rounding to the nearest integer.")
          lag <- round(lag)
        }
        lag
      },
      drift = function(drift = TRUE) {
        drift
      },
      xreg = no_xreg,
      .required_specials = c("lag")
    ),
    check = all_tsbl_checks
  )
  new_model_definition(snaive_model, !!enquo(formula), ...)
}

#' @inherit forecast.ARIMA
#' @inheritParams forecast.ETS
#' @importFrom stats qnorm time
#' @importFrom utils tail
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   forecast()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   forecast()
#' @export
forecast.RW <- function(object, new_data, specials = NULL, simulate = FALSE, bootstrap = FALSE, times = 5000, ...) {
  h <- NROW(new_data)
  lag <- object$lag
  fullperiods <- (h - 1) / lag + 1
  steps <- rep(1:fullperiods, rep(lag, fullperiods))[1:h]

  b <- object$b
  b.se <- object$b.se
  if (is_empty(b)) {
    b <- b.se <- 0
  }

  # Produce forecasts
  if (simulate || bootstrap) { # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x) {
      generate(object, new_data, bootstrap = bootstrap)[[".sim"]]
    }) %>%
      transpose() %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    distributional::dist_sample(sim)
  } else {
    fc <- rep(object$future, fullperiods)[1:h] + steps * b
    mse <- mean(residuals(object)^2, na.rm = TRUE)
    if (is.nan(mse)) mse <- NA
    # Adjust prediction intervals to allow for drift coefficient standard error
    se <- sqrt(mse * steps + (steps * b.se)^2)
    distributional::dist_normal(fc, se)
  }
}

#' @inherit generate.ETS
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   generate()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   generate()
#' @export
generate.RW <- function(x, new_data, bootstrap = FALSE, ...) {
  if (!is_regular(new_data)) {
    abort("Simulation new_data must be regularly spaced")
  }

  lag <- x$lag
  if (!is_empty(x$b)) {
    b <- stats::rnorm(1, mean = x$b, sd = x$b.se)
  } else {
    b <- 0
  }
  fits <- c(x$.fitted, x$future)

  start_idx <- min(new_data[[index_var(new_data)]])
  start_pos <- match(start_idx, seq(x$time$start, by = default_time_units(x$time$interval), length.out = length(fits)))

  future <- fits[start_pos + seq_len(lag) - 1]

  if (any(is.na(future))) {
    abort("The first lag window for simulation must be within the model's training set.")
  }

  if (!(".innov" %in% names(new_data))) {
    if (bootstrap) {
      new_data$.innov <- sample(stats::na.omit(residuals(x) - mean(residuals(x), na.rm = TRUE)),
        NROW(new_data),
        replace = TRUE
      )
    }
    else {
      new_data$.innov <- stats::rnorm(NROW(new_data), sd = sqrt(x$sigma2))
    }
  }

  sim_rw <- function(e) {
    # Cumulate errors
    dx <- e + b
    lag_grp <- rep_len(seq_len(lag), length(dx))
    dx <- split(dx, lag_grp)
    cumulative_e <- unsplit(lapply(dx, cumsum), lag_grp)
    rep_len(future, length(dx)) + cumulative_e
  }

  new_data %>%
    group_by_key() %>%
    transmute(".sim" := sim_rw(!!sym(".innov")))
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   fitted()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   fitted()
#' @export
fitted.RW <- function(object, ...) {
  object[[".fitted"]]
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   residuals()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   residuals()
#' @export
residuals.RW <- function(object, ...) {
  object[[".resid"]]
}

#' Glance a lag walk model
#'
#' Construct a single row summary of the lag walk model.
#' Contains the variance of residuals (`sigma2`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   glance()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   glance()
#' @export
glance.RW <- function(x, ...) {
  tibble(sigma2 = x[["sigma2"]])
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' as_tsibble(Nile) %>%
#'   model(NAIVE(value)) %>%
#'   tidy()
#'
#' library(tsibbledata)
#' aus_production %>%
#'   model(snaive = SNAIVE(Beer ~ lag("year"))) %>%
#'   tidy()
#' @export
tidy.RW <- function(x, ...) {
  drift <- !is_empty(x$b)
  tibble(
    term = if (drift) "b" else chr(),
    estimate = x$b, std.error = x$b.se,
    statistic = x$b / x$b.se,
    p.value = 2 * stats::pt(abs(x$b / x$b.se), length(x$.resid) - x$lag - drift, lower.tail = FALSE)
  )
}

#' @export
report.RW <- function(object, ...) {
  cat("\n")
  if (!is_empty(object[["b"]])) {
    cat(paste("Drift: ", round(object[["b"]], 4),
      " (se: ", round(object[["b.se"]], 4), ")\n",
      sep = ""
    ))
  }
  cat(paste("sigma^2:", round(object[["sigma2"]], 4), "\n"))
}

#' @importFrom stats coef
#' @export
model_sum.RW <- function(x) {
  drift <- !is_empty(x[["b"]])
  if (x[["lag"]] == 1 && !drift) {
    method <- "NAIVE"
  }
  else if (x[["lag"]] != 1) {
    method <- "SNAIVE"
  }
  else {
    method <- "RW"
  }
  if (drift) {
    method <- paste(method, "w/ drift")
  }
  method
}
