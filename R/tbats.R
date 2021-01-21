train_tbats <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by tbats.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  # Construct model output
  structure(
    list(
    ),
    class = "TBATS"
  )
}

specials_tbats <- new_specials(
  trend = function(include = TRUE, damped = FALSE) {
    as.list(environment())
  },
  transform = function(include = TRUE, lower = 0, upper = 1) {
    as.list(environment())
  },
  fourier = function(period = NULL, K = NULL) {
    as.list(environment())
  },
  season = function(period = NULL, K = NULL) {
    abort("BATS models are not yet supported.")
  },
  ARMA = function(p = 0:5, q = 0:5) {
    as.list(environment())
  },
  xreg = no_xreg,
  .required_specials = c("trend", "transform", "fourier", "ARMA")
)

#' TBATS model (Exponential smoothing state space model with Box-Cox
#' transformation, ARMA errors, Trend and Seasonal components)
#'
#' Fits a TBATS model as described in De Livera, Hyndman & Snyder (2011). 
#'
#' @references De Livera, A.M., Hyndman, R.J., & Snyder, R. D. (2011),
#' Forecasting time series with complex seasonal patterns using exponential
#' smoothing, \emph{Journal of the American Statistical Association},
#' \bold{106}(496), 1513-1527.
#'
#' @export
TBATS <- function(formula, ...) {
  tbats_model <- new_model_class("TBATS", train_tbats, specials_tbats,
    origin = NULL, check = all_tsbl_checks
  )
  new_model_definition(tbats_model, !!enquo(formula), ...)
}

#' @inherit forecast.ETS
#'
#' @export
forecast.TBATS <- function(object, new_data, specials = NULL, ...) {
  abort("Not yet implemented.")
}

#' @inherit generate.ETS
#' @export
generate.TBATS <- function(x, new_data, specials = NULL, bootstrap = FALSE, ...) {
  abort("Not yet implemented.")
}

#' @inherit fitted.ARIMA
#'
#' @export
fitted.TBATS <- function(object, ...) {
  abort("Not yet implemented.")
}

#' @inherit residuals.ARIMA
#' @export
residuals.TBATS <- function(object, ...) {
  abort("Not yet implemented.")
}

#' Glance a tbats model
#'
#' Construct a single row summary of the tbats model.
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#' @export
glance.TBATS <- function(x, ...) {
  abort("Not yet implemented.")
}

#' @inherit tidy.ARIMA
#' @export
tidy.TBATS <- function(x, ...) {
  abort("Not yet implemented.")
}

#' @export
report.TBATS <- function(object, ...) {
  abort("Not yet implemented.")
}

#' @export
model_sum.TBATS <- function(x) {
  "TBATS"
}

#' Refit a tbats model 
#'
#' Applies a fitted tbats model to a new dataset.
#'
#' @export
refit.TBATS <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  abort("Not yet implemented.")
}
