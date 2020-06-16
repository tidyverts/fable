#' @importFrom stats sd
train_theta <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by the Theta method")
  }
  
  y <- unclass(.data)[[measured_vars(.data)]]
  n <- length(y)
  
  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  # Check seasonality
  m <- specials$season[[1]]$period
  y <- ts(y, frequency = m)
  if (m > 1 && !is.constant(y) && n > 2 * m) {
    r <- as.numeric(stats::acf(y, lag.max = m, plot = FALSE)$acf)[-1]
    stat <- sqrt((1 + 2 * sum(r[-m]^2)) / n)
    if(!(abs(r[m]) / stat > qnorm(0.95))) {
      m <- 1L
    }
  } else {
    m <- 1L
  }
  
  # Seasonal decomposition
  if (m > 1L) {
    dcmp <- stats::decompose(y, type = specials$season[[1]]$method)
    if (any(abs(dcmp$seasonal) < 1e-10)) {
      warning("Seasonal indexes equal to zero. Using non-seasonal Theta method")
    } else {
      y_sa <- if(dcmp$type == "additive") dcmp$x - dcmp$seasonal else dcmp$x / dcmp$seasonal
    }
  } else {
    y_sa <- y
  }
  
  # Find theta lines
  ses <- etsmodel(y_sa, m, "A", "N", "N", FALSE, opt.crit = "mse", nmse = 3, bounds = "both")
  alpha <- pmax(1e-10, ses$par["alpha"])
  sigma2 <- sum(ses$residuals ^ 2, na.rm = TRUE) / (n - length(ses$par))
  drift <- stats::lsfit(0:(n - 1), y_sa)$coefficients[2] / 2
  
  # Reseasonalize
  if (m > 1L) {
    ses$fitted <- if(dcmp$type == "additive") ses$fitted + dcmp$seasonal else ses$fitted * dcmp$seasonal
    ses$residuals <- y - ses$fitted
  }
  
  structure(
    list(
      fitted = as.numeric(ses$fitted),
      resid = as.numeric(ses$residuals),
      period = m,
      alpha = alpha,
      l0 = ses$par["l"],
      lT = ses$states[n+1,1],
      drift = drift,
      sigma2 = sigma2,
      dcmp = dcmp$type,
      season = if(m > 1L) dcmp$seasonal[seq(n-m+1, n)] else NULL
    ),
    class = "fable_theta"
  )
}

specials_theta <- new_specials(
  season = function(period = NULL, method = c("multiplicative", "additive")) {
    period <- get_frequencies(period, self$data, .auto = "smallest")
    method <- match.arg(method)
    list(period = period, method = method)
  },
  .required_specials = "season"
)

#' Theta method
#'
#' The theta method of Assimakopoulos and Nikolopoulos (2000) is equivalent to
#' simple exponential smoothing with drift. This is demonstrated in Hyndman and
#' Billah (2003).
#'
#' The series is tested for seasonality using the test outlined in A&N. If
#' deemed seasonal, the series is seasonally adjusted using a classical
#' multiplicative decomposition before applying the theta method. The resulting
#' forecasts are then reseasonalized.
#'
#' More general theta methods are available in the forecTheta package.
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @section Specials:
#' 
#' \subsection{season}{
#' The `season` special is used to specify the parameters of the seasonal adjustment via classical decomposition.
#' \preformatted{
#' window(period = NULL, method = c("multiplicative", "additive"))
#' }
#'
#' \tabular{ll}{
#'   `period` \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").  \cr
#'   `method` \tab The type of classical decomposition to apply. The original Theta method always used multiplicative seasonal decomposition, and so this is the default.
#' }
#' }
#' 
#' @return A model specification.
#'
#' @references 
#' Assimakopoulos, V. and Nikolopoulos, K. (2000). The theta model:
#' a decomposition approach to forecasting. \emph{International Journal of
#' Forecasting} \bold{16}, 521-530.
#'
#' Hyndman, R.J., and Billah, B. (2003) Unmasking the Theta method.
#' \emph{International J. Forecasting}, \bold{19}, 287-290.
#' 
#' @author Rob J Hyndman, Mitchell O'Hara-Wild
#' @export
THETA <- function(formula, ...) {
  theta_model <- new_model_class("theta",
                                train = train_theta,
                                specials = specials_theta
  )
  new_model_definition(theta_model, !!enquo(formula), ...)
}

#' @importFrom fabletools forecast
#'
#' @inherit forecast.ARIMA
#' @export
forecast.fable_theta <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 5000, ...) {
  if (bootstrap) {
    abort("Bootstrapped forecasts are not yet supported for the Theta method.")
  }
  h <- NROW(new_data)
  
  n <- length(object$resid)
  alpha <- object$alpha
  drift <- object$drift
  sigma2 <- object$sigma2
  m <- object$period
  
  # Produce forecasts
  fc <- ets_fc_class1(h, object$lT, "N", "N", FALSE, m, sigma2, par = alpha)
  fc <- fc$mu + drift * (0:(h - 1) + (1 - (1 - alpha)^n) / alpha)
  
  # Re-seasonalise
  if(m > 1L){
    seas_fc <- rep(object$season, trunc(1 + h / m))[1:h]
    fc <- if(object$dcmp == "additive") fc + seas_fc else fc * seas_fc 
  }
  
  se <- sqrt(sigma2) * sqrt((0:(h - 1)) * alpha^2 + 1)
  distributional::dist_normal(fc, se)
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   fitted()
#' @export
fitted.fable_theta <- function(object, ...) {
  object$fitted
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   residuals()
#' @export
residuals.fable_theta <- function(object, ...) {
  object$resid
}

#' Glance a theta method
#'
#' Construct a single row summary of the average method model.
#'
#' Contains the variance of residuals (`sigma2`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#' @export
glance.fable_theta <- function(x, ...) {
  tibble(sigma2 = x$sigma2)
}

#' @inherit tidy.ARIMA
#'
#' @export
tidy.fable_theta <- function(x, ...) {
  tibble(
    term = c("alpha", "level", "drift"), 
    estimate = c(x$alpha, x$l0, x$drift)
  )
}

#' @export
report.fable_theta <- function(object, ...) {
  cat("\n")
  cat(paste("Alpha:", round(object$alpha, 4), "\n"))
  cat(paste("Drift:", round(object$drift, 4), "\n"))
  cat(paste("sigma^2:", round(object$sigma2, 4), "\n"))
}

#' @export
model_sum.fable_theta <- function(x) {
  paste0("THETA")
}
