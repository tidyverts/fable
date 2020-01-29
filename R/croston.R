#' Croston's method
#' 
#' Based on Croston's (1972) method for intermittent demand forecasting, also described in Shenstone and Hyndman (2005). Croston's method involves using simple exponential smoothing (SES) on the non-zero elements of the time series and a separate application of SES to the times between non-zero elements of the time series. 
#' 
#' Note that forecast distributions are not computed as Croston's method has no underlying stochastic model. In a later update, we plan to support distributions via the equivalent stochastic models that underly Croston's method (Shenstone and Hyndman, 2005)
#' 
#' @param formula Model specification (see "Specials" section).
#' @param opt_crit The optimisation criterion used to optimise the parameters.
#' @param ... Not used.
#' 
#' @section Specials:
#' 
#' \subsection{demand}{
#' The `demand` special specifies parameters for the demand SES application.
#' \preformatted{
#' demand(initial = NULL, param = NULL)
#' }
#'
#' \tabular{ll}{
#'   `initial`  \tab The initial value for the demand application of SES. \cr
#'   `param`    \tab The smoothing parameter for the demand application of SES.
#' }
#' }
#' 
#' \subsection{interval}{
#' The `interval` special specifies parameters for the interval SES application.
#' \preformatted{
#' interval(initial = NULL, param = NULL)
#' }
#'
#' \tabular{ll}{
#'   `initial`  \tab The initial value for the interval application of SES. \cr
#'   `param`    \tab The smoothing parameter for the interval application of SES.
#' }
#' }
#' 
#' @return A model specification.
#' 
#' @references Croston, J. (1972) "Forecasting and stock control for
#' intermittent demands", \emph{Operational Research Quarterly}, \bold{23}(3),
#' 289-303.
#'
#' Shenstone, L., and Hyndman, R.J. (2005) "Stochastic models underlying
#' Croston's method for intermittent demand forecasting". \emph{Journal of
#' Forecasting}, \bold{24}, 389-402.
#' 
#' Kourentzes, N. (2014) "On intermittent demand model optimisation and 
#' selection". \emph{International Journal of Production Economics}, \bold{156},
#' 180-190. \url{http://dx.doi.org/10.1016/j.ijpe.2014.06.007}.
#' 
#' @examples 
#' library(tsibble)
#' poisson <- tsibble(
#'   time = yearmonth("2012 Dec") + seq_len(24),
#'   count = rpois(24, lambda = 0.3)
#' )
#' 
#' poisson %>% 
#'   autoplot(count)
#' 
#' poisson %>% 
#'   model(CROSTON(count)) %>% 
#'   forecast(h = "2 years") %>% 
#'   autoplot(poisson)
#' 
#' @importFrom stats model.matrix
#' @export
CROSTON <- function(formula, opt_crit = c("mse", "mae"), ...){
  opt_crit <- match.arg(opt_crit)
  croston_model <- new_model_class("croston", train = train_croston, 
                                   specials = specials_croston,
                                   check = all_tsbl_checks)
  new_model_definition(croston_model, !!enquo(formula), opt_crit = opt_crit, ...)
}


specials_croston <- new_specials(
  demand = function(initial = NULL, param = NULL){
    
    if(!is.null(initial) && initial < 0){
      abort("The initial demand for Croston's method must be non-negative")
    }
    
    as.list(environment())
  },
  interval = function(initial = NULL, param = NULL, method = c("mean", "naive")){
    method <- match.arg(method)
    
    if(!is.null(initial) && initial < 1){
      abort("The initial interval for Croston's method must be greater than (or equal to) 1.")
    }
    
    as.list(environment())
  },
  .required_specials = c("demand", "interval")
)

train_croston <- function(.data, specials, opt_crit = "mse", ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by Croston's method.")
  }
  
  # Get response
  y <- unclass(.data)[[measured_vars(.data)]]
  
  # Check data
  if(any(y < 0)){
    abort("All observations must be non-negative for Croston's method.")
  }
  
  non_zero <- which(y != 0)
  
  if(length(non_zero)<2){
    abort("At least two non-zero values are required to use Croston's method.")
  }
  
  # Get specials
  demand <- specials$demand[[1]]
  interval <- specials$interval[[1]]
  
  # Croston demand/interval decomposition
  y_demand <- y[non_zero]
  y_interval <- c(non_zero[1], diff(non_zero))
  
  # Initialise parameters
  if(is.null(demand$initial)){
    demand$initial <- y_demand[1]
  }
  if(is.null(interval$initial)){
    interval$initial <- switch(interval$method, mean = mean(y_interval), naive = y_interval[1])
  }
  
  # Optimise parameters
  # Smoothing params init at 0.05
  # par <- c(demand$initial, interval$initial, demand$param, interval$param)
  par <- c(demand$initial, interval$initial, 0.05, 0.05)
  par <- stats::optim(
    par=par, optim_croston, demand = demand, interval = interval,
    y = y, y_demand = y_demand, y_interval = y_interval, 
    non_zero = non_zero, n = length(y), opt_crit = opt_crit,
    lower = c(0, 1, 0, 0), upper = c(max(y_demand), max(y_interval), 1, 1),
    method="L-BFGS-B", control = list(maxit=2000)
  )$par
  
  demand$initial <- par[1]
  interval$initial <- par[2]
  demand$param <- par[3]
  interval$param <- par[4]
  
  fit <- estimate_croston(y_demand, y_interval, demand, interval, 
                          non_zero = non_zero, n = length(y))
  
  structure(
    list(
      par = list(
        term = c("demand_init", "demand_par", "interval_init", "interval_par"),
        estimate = c(demand$initial, demand$param, interval$initial, interval$param)
      ),
      .fitted = fit,
      .resid = y - fit
    ),
    class = "croston"
  )
}

optim_croston <- function(par, demand, interval, y, opt_crit, ...){
  demand$initial <- par[1]
  interval$initial <- par[2]
  demand$param <- par[3]
  interval$param <- par[4]
  
  frc.in <- estimate_croston(..., demand = demand, interval = interval)
  resid <- y - frc.in
  
  if (opt_crit == "mse"){
    mean(resid^2, na.rm = TRUE)
  } else if(opt_crit == "mae"){
    mean(abs(resid), na.rm = TRUE)
  }
}

estimate_croston <- function(y_demand, y_interval, demand, interval, non_zero, n){
  # Initialise fits
  k <- length(y_demand)
  fit_demand <- numeric(k)
  fit_interval <- numeric(k)
  fit_demand[1] <- demand$initial
  fit_interval[1] <- interval$initial
  
  # Fit model
  for (i in 2:k){
    fit_demand[i] <- fit_demand[i-1] + demand$param * (y_demand[i] - fit_demand[i-1])
    fit_interval[i] <- fit_interval[i-1] + interval$param * (y_interval[i] - fit_interval[i-1])
  }
  ratio <- fit_demand/fit_interval
  
  # In-sample demand rate
  rep(c(NA_real_, ratio), diff(c(0, non_zero, n)))
}

#' @inherit forecast.ARIMA
#'
#' @examples 
#' library(tsibble)
#' poisson <- tsibble(
#'   time = yearmonth("2012 Dec") + seq_len(24),
#'   count = rpois(24, lambda = 0.3)
#' )
#' 
#' poisson %>% 
#'   model(CROSTON(count)) %>% 
#'   forecast()
#'
#' @export
forecast.croston <- function(object, new_data, specials = NULL, ...){
  h <- nrow(new_data)
  fc <- rep(object$.fitted[length(object$.fitted)], h)
  construct_fc(fc, numeric(h), dist_unknown(h))
}

#' @inherit fitted.ARIMA
#' 
#' @examples 
#' library(tsibble)
#' poisson <- tsibble(
#'   time = yearmonth("2012 Dec") + seq_len(24),
#'   count = rpois(24, lambda = 0.3)
#' )
#' 
#' poisson %>% 
#'   model(CROSTON(count)) %>% 
#'   tidy()
#' @export
fitted.croston <- function(object, ...){
  object[[".fitted"]]
}


#' @inherit residuals.ARIMA
#' 
#' @examples 
#' library(tsibble)
#' poisson <- tsibble(
#'   time = yearmonth("2012 Dec") + seq_len(24),
#'   count = rpois(24, lambda = 0.3)
#' )
#' 
#' poisson %>% 
#'   model(CROSTON(count)) %>% 
#'   residuals()
#' @export
residuals.croston <- function(object, ...){
  object[[".resid"]]
}

#' @inherit tidy.ARIMA
#' 
#' @examples 
#' library(tsibble)
#' poisson <- tsibble(
#'   time = yearmonth("2012 Dec") + seq_len(24),
#'   count = rpois(24, lambda = 0.3)
#' )
#' 
#' poisson %>% 
#'   model(CROSTON(count)) %>% 
#'   tidy()
#'   
#' @export
#' @export
tidy.croston <- function(x, ...){
  as_tibble(x[["par"]])
}