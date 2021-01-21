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

fit_tbats <- function(y, m, trend, dampen, transform, ar, ma) {
  # Test inputs
  # y = USAccDeaths
  # m = 12
  # k = 3
  # trend = TRUE
  # dampen = TRUE
  # transform = TRUE
  # ar = ma = numeric()
  
  
  p = length(ar)
  q = length(ma)
  
  par <- initial_tbats_par(y, m, k, trend, dampen, transform)
  
  # Initialise x0 matrix as x0=0 for l, b, s (2*k), p, q
  x0 <- matrix(c(0, numeric(as.integer(is.null(par$beta))), numeric(2*sum(k)), numeric(p), numeric(q)))
  
  # Create transposed w matrix
  # Inputs represent: level, dampening, fourier harmonics, ar, ma
  wt <- matrix(c(1, par$phi, rep(rep(c(1,0), length(k)), rep(k, each = 2)), ar, ma),
               nrow = 1)
  
  # Create gamma vector
  # gamma = c(rep(gamma1[1], k[1]), rep(gamma2[1], k[1]), rep(gamma1[2], k[2]), rep(gamma2[2], k[2]), ...)
  gamma <- matrix(rep(cbind(par$gamma1,par$gamma2), rep(k, length(k))),
                  nrow = 1)
  
  # Create g vector
  # Inputs represent: alpha, beta, gamma, <AR(p): 1, 0*p-1>, <MA(q): 1, 0*q-1>
  g <- matrix(
    c(par$alpha, par$beta, gamma, 
      1[p>0], numeric(max(0, p-1)),
      1[q>0], numeric(max(0, q-1))),
    ncol = 1)
  
  # Create F matrix
  Fmat <- matrix(0, nrow = nrow(g), ncol = nrow(g))
  Fmat[1,1] <- 1
  # to be continued...
}

initial_tbats_par <- function(y, m, k, trend, dampen, transform) {
  alpha <- 0.09
  beta <- gamma1 <- gamma2 <- lambda <- NULL
  if (trend) {
    beta <- 0.05
    if (dampen) {
      phi <- .999
    } else {
      phi <- 1
    }
  }
  if (!is_empty(m)) {
    gamma1 <- rep(0, length(k))
    gamma2 <- rep(0, length(k))
  }
  if (transform) {
    require_package("feasts")
    lambda <- unname(feasts::guerrero(y, lower = 0, upper = 1.5, .period = m))
  }
  list(alpha = alpha, beta = beta, phi = phi,
       gamma1 = gamma1, gamma2 = gamma2, lambda = lambda)
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
