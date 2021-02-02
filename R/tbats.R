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
  y = USAccDeaths
  m = 12
  k = 3
  tau <- sum(k)*2
  trend = TRUE
  dampen = TRUE
  transform = TRUE
  bc_lambda = NULL
  ar = ma = numeric()
  
  p = length(ar)
  q = length(ma)
  
  par <- initial_tbats_par(y, m, k, trend, dampen, transform)
  
  # Initialise x0 matrix as x0=0 for l, b, s (2*k), p, q
  x0 <- matrix(c(0, numeric(as.integer(!is.null(par$beta))), numeric(tau), numeric(p), numeric(q)))
  
  # Create transposed w matrix
  # Inputs represent: level, dampening, fourier harmonics, ar, ma
  wt <- matrix(c(1, par$phi, rep(rep(c(1,0), length(k)), rep(k, each = 2)), ar, ma),
               nrow = 1)
  
  # Create gamma vector
  # gamma = c(rep(gamma1[1], k[1]), rep(gamma2[1], k[1]), rep(gamma1[2], k[2]), rep(gamma2[2], k[2]), ...)
  gamma <- rep(cbind(par$gamma1,par$gamma2), rep(k, length(k)))
  if(!is.null(gamma)) gamma <- matrix(gamma, nrow = 1)
  
  # Create g vector
  # Inputs represent: alpha, beta, gamma, <AR(p): 1, 0*p-1>, <MA(q): 1, 0*q-1>
  g <- matrix(
    c(par$alpha, par$beta, gamma, 
      c(1)[p>0], numeric(max(0, p-1)),
      c(1)[q>0], numeric(max(0, q-1))
    ),
    ncol = 1)
  
  Fmat <- make_tbats_Fmat(par$alpha, par$beta, par$phi, m, k, gamma, ar, ma)
  
  Dmat <- F - g%*%wt
  
  y_trans <- if(is.null(bc_lambda)) y else box_cox(y, lambda = bc_lambda) 
  
  # Calculate (initial?) yhat - \tilde{y}
  # y_tilde <- filter_tbats(y_trans, x0, Fmat, g, wt)
  
  ## Initialise states (and 'find the seed states'?)
  len <- length(y_trans)
  x <- matrix(nrow = nrow(x0), ncol = len)
  y_hat <- y_tilde <- matrix(nrow = 1, ncol = len)
  y_hat[1] <- wt %*% x0
  y_tilde[1] <- y_trans[1] - y_hat[1]
  x[,1] <- Fmat %*% x0 + g%*%y_tilde[1]
  ## Filter series
  for (t in 2:len) {
    y_hat[t] <- wt %*% x[,t-1]
    y_tilde[t] <- y_trans[t] - y_hat[t]
    x[,t] <- Fmat %*% x[,t-1] + g%*%y_tilde[t]
  }
  
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

make_tbats_Fmat <- function(alpha, beta, phi, m, k, gamma, ar, ma) {
  p <- length(ar)
  q <- length(ma)
  tau <- length(gamma)
  # Create F matrix
  Fmat <- matrix(0, nrow = nrow(g), ncol = nrow(g))
  # level
  Fmat[1,1] <- 1
  # trend
  if(!is.null(par$beta)) {
    Fmat[c(2, 1:2 + nrow(Fmat))] <- c(0, par$phi, par$phi)
  }
  # fourier
  if(!is.null(gamma)) {
    Amat <- matrix(0, tau, tau)
    fourier_mat <- function(m, k) {
      Ci <- diag(cos(2*pi*seq_len(k)/m))
      Si <- diag(sin(2*pi*seq_len(k)/m))
      Ai <- matrix(nrow = k*2, ncol = k*2)
      pos1 <- seq_len(k)
      pos2 <- seq(k+1, length.out = k)
      Ai[pos1, pos1] <- Ci
      Ai[pos1, pos2] <- Si
      Ai[pos2, pos1] <- -Si
      Ai[pos2, pos2] <- Ci
      Ai
    }
    idx <- c(1,cumsum(2*k))
    for(i in seq_along(m)) {
      Amat[idx[i]:idx[i+1],idx[i]:idx[i+1]] <- fourier_mat(m[i], k[i])
    }
    Aloc <- seq(2 + !is.null(par$beta), length.out = tau)
    Fmat[Aloc, Aloc] <-Amat
  }
  # ar
  if(p > 0) {
    ar_cols <- seq(2 + tau + !is.null(par$beta), length.out = p)
    
    # Add ar:level
    Fmat[1, ar_cols] <- par$alpha * ar # alpha * phi
    cur_row <- 2
    # Add ar:trend
    if(!is.null(par$beta)) {
      Fmat[cur_row, ar_cols] <- par$beta * ar # beta * phi
      cur_row <- cur_row + 1
    }
    # Add ar:fourier
    if(!is.null(gamma)) {
      Fmat[seq(cur_row, length.out = tau), ar_cols] <- t(gamma) %*% ar # B matrix
      cur_row <- cur_row + tau
    }
    # Add ar component
    Fmat[cur_row, ar_cols] <- ar
    Fmat[seq(cur_row + 1, length.out = p-1), ar_cols] <- diag(nrow = p-1, ncol = p) # I<p-1,p>
  }
  # ma
  if(q > 0) {
    ma_cols <- seq(2 + tau + p + !is.null(par$beta), length.out = q)
    
    # Add ma:level
    Fmat[1, ma_cols] <- par$alpha * ma # alpha * theta
    cur_row <- 2
    # Add ma:trend
    if(!is.null(par$beta)) {
      Fmat[cur_row, ma_cols] <- par$beta * ma # beta * theta
      cur_row <- cur_row + 1
    }
    # Add ar:fourier
    if(!is.null(gamma)) {
      Fmat[seq(cur_row, length.out = tau), ma_cols] <- t(gamma) %*% ma # C matrix
      cur_row <- cur_row + tau
    }
    # Add ar component
    Fmat[cur_row, ma_cols] <- ma
    Fmat[seq(cur_row + p + 1, length.out = q-1), ma_cols] <- diag(nrow = q-1, ncol = q) # I<p-1,p>
  }
  Fmat
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
