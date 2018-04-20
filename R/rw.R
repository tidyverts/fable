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
#' @param data A data frame
#' @param formula Model specification.
#' 
#' @examples 
#' elecdemand %>% RW(Demand ~ drift())
#' 
#' @export
RW <- function(data, formula = ~ lag(1)){
  # Capture call
  cl <- new_quosure(match.call())
  formula <- enquo(formula)
  
  # Coerce data
  data <- as_tsibble(data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    lag = function(lag = 1){
      list(lag = lag)
    },
    drift = function(drift = TRUE){
      list(include.constant = drift)
    },
    xreg = specials_xreg,
    
    parent_env = get_env(cl)
  )
  
  # Parse model
  model <- data %>% 
    parse_model(formula, specials = specials)
  
  # Output model
  eval_tidy(expr(model_RW(data, model, !!!flatten_first_args(model$args))))
}

model_RW <- function(data, model, lag = 1, ...){
  lag <- get_frequencies(lag, data)
  
  lag_arg <- if(lag == 1){list(order = c(0, 1, 0))} else {list(seasonal = list(order = c(0, 1, 0), period = lag))}
  out <- wrap_ts_model(data, "Arima", model, !!!lag_arg, period = lag, ...)
  out[["model"]][[1]] <- enclass(out[["model"]][[1]], "RW")
  out
  # TODO: Adjust prediction intervals to allow for drift coefficient standard error
}

#' @rdname RW
#'
#' @examples
#' 
#' forecast::gold %>% NAIVE
#'
#' @export
NAIVE <- RW

#' @rdname RW
#'
#' @examples
#' 
#' elecdemand %>% SNAIVE(Temperature ~ lag("day"))
#'
#' @export
SNAIVE <- function(data, formula = ~ lag("smallest")){
  # Capture call
  cl <- new_quosure(match.call())
  formula <- enquo(formula)
  
  # Coerce data
  data <- as_tsibble(data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    lag = function(lag = 1){
      list(lag = lag)
    },
    drift = function(drift = TRUE){
      list(include.constant = drift)
    },
    xreg = specials_xreg,
    
    parent_env = get_env(cl)
  )
  
  # Parse model
  model <- data %>% 
    parse_model(formula, specials = specials)
  
  # Output model
  eval_tidy(expr(model_SNAIVE(data, model, !!!flatten_first_args(model$args))))
}

model_SNAIVE <- function(data, model, lag = "smallest", ...){
  lag <- get_frequencies(lag, data)
  if(lag == 1){
    abort("Non-seasonal model specification provided, use RW() or provide a different lag specification.")
  }
  
  model_RW(data, model, lag, ...)
}

#' @importFrom dplyr case_when
model_sum.RW <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  drift <- NCOL(x$xreg) == 1 && is.element("drift", names(x$coef))
  result <- case_when(
    order[2]==1 && drift ~ "RW",
    order[2]==1 ~ "NAIVE",
    order[5]==1 ~ "SNAIVE",
    TRUE ~ "RW"
    )
  
  if (!is.null(x$xreg)) {
    if (drift) {
      result <- paste(result, "w/ drift")
    } else {
      result <- paste("Regression with", result, "errors")
    }
  }
  result
}