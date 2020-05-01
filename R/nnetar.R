#' @importFrom stats ar complete.cases
train_nnetar <- function(.data, specials, n_nodes, n_networks, scale_inputs, ...) {
  require_package("nnet")

  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by NNETAR.")
  }

  y <- x <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }

  n <- length(x)

  if (n < 3) {
    stop("Not enough data to fit a model")
  }

  # Get args
  p <- P <- period <- NULL
  assignSpecials(specials["AR"])
  xreg <- specials$xreg[[1]]

  # Check for constant data in time series
  constant_data <- is.constant(x)
  if (constant_data) {
    warn("Constant data, setting `AR(p=1, P=0)`, and `scale_inputs=FALSE`")
    scale_inputs <- FALSE
    p <- 1
    P <- 0
  }
  ## Check for constant data in xreg
  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (any(apply(xreg, 2, is.constant))) {
      warn("Constant xreg column, setting `scale_inputs=FALSE`")
      scale_inputs <- FALSE
    }
  }

  # Scale inputs
  y_scale <- xreg_scale <- NULL
  if (scale_inputs) {
    x <- scale(x, center = TRUE, scale = TRUE)
    # y_scale <- tibble(
    #   term = c("y_center", "y_scale"),
    #   estimate = c(attr(x, "scaled:center"), attr(x, "scaled:scale"))
    # )
    y_scale <- list(
      center = attr(x, "scaled:center"),
      scale = attr(x, "scaled:scale")
    )
    x <- c(x)

    if (!is.null(xreg)) {
      xreg <- scale(xreg, center = TRUE, scale = TRUE)
      # xreg_scale <- tibble(
      #   term = c("xreg_center", "xreg_scale"),
      #   estimate = c(attr(xreg, "scaled:center"), attr(xreg, "scaled:scale"))
      # )
      xreg_scale <- list(
        center = attr(xreg, "scaled:center"),
        scale = attr(xreg, "scaled:scale")
      )
    }
  }

  # Construct lagged matrix
  if (is.null(p)) {
    if (period > 1 && n > 2 * period) {
      y_sa <- stats::stl(ts(x, frequency = period), s.window = 13)
      y_sa <- y_sa$time.series[, "trend"] + y_sa$time.series[, "remainder"]
    }
    else {
      y_sa <- x
    }
    p <- max(length(ar(y_sa)$ar), 1)
  }
  if (p >= n) {
    warn("Reducing number of lagged inputs due to short series")
    p <- n - 1
  }
  lags <- 1:p
  if (P > 0 && n >= period * P + 2) {
    lags <- sort(unique(c(lags, period * (1:P))))
  }

  maxlag <- max(lags)
  nlag <- length(lags)
  x_lags <- matrix(NA_real_, ncol = nlag, nrow = n - maxlag)
  for (i in 1:nlag) {
    x_lags[, i] <- x[(maxlag - lags[i] + 1):(n - lags[i])]
  }
  x <- x[-(1:maxlag)]
  # Add xreg into lagged matrix
  x_lags <- cbind(x_lags, xreg[-(1:maxlag), ])
  if (is.null(n_nodes)) {
    n_nodes <- round((NCOL(x_lags) + 1) / 2)
  }

  # Remove missing values if present
  j <- complete.cases(x_lags, x)
  x_lags <- x_lags[j, , drop = FALSE]
  x <- x[j]
  ## Stop if there's no data to fit
  if (NROW(x_lags) == 0) {
    abort("No data to fit (possibly due to missing values)")
  }

  # Fit the nnet
  nn_models <- map(
    seq_len(n_networks),
    function(.) wrap_nnet(x_lags, x, size = n_nodes, ...)
  )

  # Calculate fitted values
  pred <- map(nn_models, predict) %>%
    transpose() %>%
    map_dbl(function(x) mean(unlist(x)))
  fits <- rep(NA_real_, length(y))
  fits_idx <- c(rep(FALSE, maxlag), j)
  fits[fits_idx] <- pred
  if (scale_inputs) {
    fits <- fits * y_scale$scale + y_scale$center
  }
  res <- y - fits

  # Construct model output
  structure(
    list(
      model = nn_models,
      par = tibble(),
      est = tibble(.fitted = fits, .resid = res),
      fit = tibble(sigma2 = stats::var(res, na.rm = TRUE)),
      spec = tibble(period = period, p = p, P = P, size = n_nodes, lags = list(lags)),
      scales = list(y = y_scale, xreg = xreg_scale),
      future = mutate(
        new_data(.data, maxlag),
        !!expr_text(model_lhs(self)) := utils::tail(x, maxlag)
      )
    ),
    class = "NNETAR"
  )
}

# Wrap nnet to change the default for linout to be TRUE
wrap_nnet <- function(x, y, linout = TRUE, trace = FALSE, ...) {
  nnet::nnet(x, y, linout = linout, trace = trace, ...)
}

specials_nnetar <- new_specials(
  AR = function(p = NULL, P = 1, period = NULL) {
    period <- get_frequencies(period, self$data, .auto = "smallest")
    if (period == 1) {
      if (!missing(P) && P > 0) {
        warn("Non-seasonal data, ignoring seasonal lags")
      }
      P <- 0
    }
    if (P > 0 && NROW(self$data) < period * P + 2) {
      warn("Series too short for seasonal lags")
      P <- 0
    }
    as.list(environment())
  },
  common_xregs,
  xreg = model_xreg,
  .required_specials = c("AR"),
  .xreg_specials = names(common_xregs)
)

#' Neural Network Time Series Forecasts
#'
#' Feed-forward neural networks with a single hidden layer and lagged inputs
#' for forecasting univariate time series.
#'
#' A feed-forward neural network is fitted with lagged values of the response as
#' inputs and a single hidden layer with `size` nodes. The inputs are for
#' lags 1 to `p`, and lags `m` to `mP` where
#' `m` is the seasonal period specified.
#'
#' If exogenous regressors are provided, its columns are also used as inputs.
#' Missing values are currently not supported by this model.
#' A total of `repeats` networks are
#' fitted, each with random starting weights. These are then averaged when
#' computing forecasts. The network is trained for one-step forecasting.
#' Multi-step forecasts are computed recursively.
#'
#' For non-seasonal data, the fitted model is denoted as an NNAR(p,k) model,
#' where k is the number of hidden nodes. This is analogous to an AR(p) model
#' but with non-linear functions. For seasonal data, the fitted model is called
#' an NNAR(p,P,k)\[m\] model, which is analogous to an ARIMA(p,0,0)(P,0,0)\[m\]
#' model but with non-linear functions.
#'
#' @aliases report.NNETAR
#'
#' @param formula Model specification (see "Specials" section).
#' @param n_nodes Number of nodes in the hidden layer. Default is half of the
#' number of input nodes (including external regressors, if given) plus 1.
#' @param n_networks Number of networks to fit with different random starting
#' weights. These are then averaged when producing forecasts.
#' @param scale_inputs If TRUE, inputs are scaled by subtracting the column
#' means and dividing by their respective standard deviations. Scaling is
#' applied after transformations.
#' @param ... Other arguments passed to `\link[nnet]{nnet}`.
#'
#' @section Specials:
#'
#' \subsection{AR}{
#' The `AR` special is used to specify auto-regressive components in each of the
#' nodes of the neural network.
#'
#' \preformatted{
#' AR(p = NULL, P = 1, period = NULL)
#' }
#'
#' \tabular{ll}{
#'   `p`        \tab The order of the non-seasonal auto-regressive (AR) terms. If `p = NULL`, an optimal number of lags will be selected for a linear AR(p) model via AIC. For seasonal time series, this will be computed on the seasonally adjusted data (via STL decomposition). \cr
#'   `P`        \tab The order of the seasonal auto-regressive (SAR) terms. \cr
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").
#' }
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an NNETAR model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#' \preformatted{
#' xreg(...)
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)
#' }
#' }
#'
#' @return A model specification.
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15)))
#' @seealso
#' [Forecasting: Principles and Practices, Neural network models (section 11.3)](https://otexts.com/fpp2/nnetar.html)
#'
#' @export
NNETAR <- function(formula, n_nodes = NULL, n_networks = 20, scale_inputs = TRUE, ...) {
  nnetar_model <- new_model_class("NNETAR", train_nnetar, specials_nnetar,
    origin = NULL, check = all_tsbl_checks
  )
  new_model_definition(nnetar_model, !!enquo(formula),
    n_nodes = n_nodes,
    n_networks = n_networks, scale_inputs = scale_inputs, ...
  )
}

#' @inherit forecast.ETS
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   forecast(times = 10)
#' @export
forecast.NNETAR <- function(object, new_data, specials = NULL, simulate = TRUE, bootstrap = FALSE, times = 1000, ...) {
  require_package("nnet")

  # Prepare xreg
  xreg <- specials$xreg[[1]]

  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (!is.null(object$scales$xreg)) {
      xreg <- scale(xreg, center = object$scales$xreg$center, scale = object$scales$xreg$scale)
    }
  }

  # Extract model attributes
  lags <- object$spec$lags[[1]]
  maxlag <- max(lags)
  future_lags <- rev(object$future[[measured_vars(object$future)]])

  # Compute forecast intervals
  if (!simulate) {
    warn("Analytical forecast distributions are not available for NNETAR.")
    times <- 0
  }
  sim <- map(seq_len(times), function(x) {
    generate(object, new_data, specials = specials, bootstrap = bootstrap)[[".sim"]]
  })
  if (length(sim) > 0) {
    sim <- sim %>%
      transpose() %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    distributional::dist_sample(sim)
  }
  else {
    # Compute 1-step forecasts iteratively
    h <- NROW(new_data)
    fc <- numeric(h)
    for (i in seq_len(h))
    {
      fcdata <- c(future_lags[lags], xreg[i, ])
      if (any(is.na(fcdata))) {
        abort("I can't use NNETAR to forecast with missing values near the end of the series.")
      }
      fc[i] <- mean(map_dbl(object$model, predict, newdata = fcdata))
      future_lags <- c(fc[i], future_lags[-maxlag])
    }
    # Re-scale point forecasts
    if (!is.null(object$scales$y)) {
      fc <- fc * object$scales$y$scale + object$scales$y$center
    }
    distributional::dist_degenerate(fc)
  }
}

#' @inherit generate.ETS
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   generate()
#' @export
generate.NNETAR <- function(x, new_data, specials = NULL, bootstrap = FALSE, ...) {
  # Prepare xreg
  xreg <- specials$xreg[[1]]

  if (!is.null(xreg)) {
    xreg <- as.matrix(xreg)
    if (!is.null(x$scales$xreg)) {
      xreg <- scale(xreg, center = x$scales$xreg$center, scale = x$scales$xreg$scale)
    }
  }

  if (!(".innov" %in% names(new_data))) {
    if (bootstrap) {
      res <- stats::na.omit(x$est[[".resid"]] - mean(x$est[[".resid"]], na.rm = TRUE))
      if (!is.null(x$scales$y)) {
        res <- res / x$scales$y$scale
      }
      new_data$.innov <- sample(res, NROW(new_data), replace = TRUE)
    }
    else {
      if (!is.null(x$scales$y)) {
        sigma <- sd(x$est[[".resid"]] / x$scales$y$scale, na.rm = TRUE)
      }
      else {
        sigma <- x$fit$sigma
      }
      new_data$.innov <- stats::rnorm(NROW(new_data), sd = sigma)
    }
  }
  else {
    if (!is.null(x$scales$y)) {
      new_data[[".innov"]] <- new_data[[".innov"]] / x$scales$y$scale
    }
  }

  # Extract model attributes
  lags <- x$spec$lags[[1]]
  maxlag <- max(lags)
  future_lags <- rev(x$future[[measured_vars(x$future)]])

  sim_nnetar <- function(e) {
    path <- numeric(length(e))
    for (i in 1:length(e))
    {
      fcdata <- c(future_lags[lags], xreg[i, ])
      if (any(is.na(fcdata))) {
        abort("I can't use NNETAR to forecast with missing values near the end of the series.")
      }
      path[i] <- mean(map_dbl(x$model, predict, newdata = fcdata)) + e[i]
      future_lags <- c(path[i], future_lags[-maxlag])
    }

    # Re-scale simulated paths
    if (!is.null(x$scales$y)) {
      path <- path * x$scales$y$scale + x$scales$y$center
    }
    path
  }

  new_data %>%
    group_by_key() %>%
    transmute(".sim" := sim_nnetar(!!sym(".innov")))
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   fitted()
#' @export
fitted.NNETAR <- function(object, ...) {
  object$est[[".fitted"]]
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   residuals()
#' @export
residuals.NNETAR <- function(object, ...) {
  object$est[[".resid"]]
}

#' Glance a NNETAR model
#'
#' Construct a single row summary of the NNETAR model.
#' Contains the variance of residuals (`sigma2`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   glance()
#' @export
glance.NNETAR <- function(x, ...) {
  x$fit
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' as_tsibble(airmiles) %>%
#'   model(nn = NNETAR(box_cox(value, 0.15))) %>%
#'   tidy()
#' @export
tidy.NNETAR <- function(x, ...) {
  x$par
}

#' @export
report.NNETAR <- function(object, ...) {
  cat(paste("\nAverage of", length(object$model), "networks, each of which is\n"))
  print(object$model[[1]])
  cat(
    "\nsigma^2 estimated as ",
    format(mean(residuals(object)^2, na.rm = TRUE), digits = 4),
    "\n",
    sep = ""
  )
  invisible(object)
}

#' @importFrom stats coef
#' @export
model_sum.NNETAR <- function(x) {
  p <- x$spec$p
  P <- x$spec$P
  sprintf(
    "NNAR(%s,%i)%s",
    ifelse(P > 0, paste0(p, ",", P), p),
    x$spec$size,
    ifelse(P > 0, paste0("[", x$spec$period, "]"), "")
  )
}
