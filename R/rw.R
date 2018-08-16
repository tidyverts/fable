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
#' library(tsibbledata)
#' elecdemand %>% RW(Demand ~ drift())
#' 
#' @export
RW <- function(data, formula = ~ lag(1)){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Coerce data
  data <- as_tsibble(data)
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    lag = function(lag = 1){
      lag <- get_frequencies(lag, .data)
      if(lag == 1){
        list(order = c(0L, 1L, 0L), seasonal = list(order = c(0L, 0L, 0L), period = 1))
      } else {
        list(order = c(0L, 0L, 0L), seasonal = list(order = c(0L, 1L, 0L), period = lag))
      }
    },
    drift = function(drift = TRUE){
      as.matrix(`colnames<-`(trend(.data, origin = origin), "drift"))
    },
    xreg = function(...){
      as_model_matrix(tibble(...))
    },
    .env = caller_env(),
    .required_specials = c("lag"),
    .vals = list(.data = data, origin = min(data[[expr_text(index(data))]]))
  )
  
  estimate_RW(data = data, formula = formula, specials = specials, cl = cl)
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
#' library(tsibbledata)
#' elecdemand %>% SNAIVE(Temperature ~ lag("day"))
#'
#' @export
SNAIVE <- function(data, formula = ~ lag("smallest")){
  # Capture user call
  cl <- call_standardise(match.call())
  # Coerce data
  data <- as_tsibble(data)
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    lag = function(lag = "smallest"){
      lag <- get_frequencies(lag, .data)
      if(lag == 1){
        abort("Non-seasonal model specification provided, use RW() or provide a different lag specification.")
      } else {
        list(order = c(0L, 0L, 0L), seasonal = list(order = c(0L, 1L, 0L), period = lag))
      }
    },
    drift = function(drift = TRUE){
      as.matrix(`colnames<-`(trend(.data, origin = origin), "drift"))
    },
    xreg = function(...){
      as_model_matrix(tibble(...))
    },
    
    .env = caller_env(),
    .required_specials = c("lag"),
    .vals = list(.data = data)
  )
  
  estimate_RW(data=data, formula = formula, specials = specials, cl = cl)
}

estimate_RW <- function(data, formula, specials, cl){
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  y <- eval_tidy(model_lhs(model_inputs$model), data = data)
  
  args <- model_inputs$args
  xreg <- cbind(args$drift[[1]], args$xreg[[1]])
  fit <- arima(
    y,
    order = args$lag[[1]]$order,
    seasonal = args$lag[[1]]$seasonal,
    xreg = xreg
  )
  
  fit$residuals[seq_len(args$lag[[1]]$seasonal$period)] <- NA
  fit$fitted <- y - fit$residuals
  
  fit$call <- cl
  
  # Output model
  mable(
    data,
    model = add_class(fit, "RW"),
    model_inputs
  )
}


#' @importFrom forecast forecast
#' @importFrom purrr map2
#' @importFrom stats qnorm time 
#' @importFrom utils tail
#' @importFrom dplyr pull
#' @export
forecast.RW <- function(object, data, h = NULL, newdata = NULL, ...){
  if(is.null(newdata)){
    if(is.null(h)){
      h <- get_frequencies("all", data) %>%
        .[.>2] %>%
        min
    }
    future_idx <- data %>% pull(!!index(.)) %>% fc_idx(h)
    newdata <- tsibble(!!!set_names(list(future_idx), expr_text(index(data))), index = !!index(data))
  }
  
  if("drift" %in% names(coef(object))){
    drift <- as.matrix(`colnames<-`(trend(newdata, origin = min(data[[expr_text(index(data))]])), "drift"))
  }
  else{
    drift <- NULL
  }
  
  xreg <- names(coef(object))[-match("drift", names(coef(object)))]
  if(length(xreg) > 0){
    xreg <- as_model_matrix(eval_tidy(tibble(!!!parse_exprs(xreg)), data = newdata))
  }
  else{
    xreg <- NULL
  }
  
  xreg <- cbind(drift, xreg)
  
  object$call$xreg <- xreg # Bypass predict.Arima NCOL check
  fc <- predict(object, n.ahead = NROW(newdata), newxreg = xreg, ...)
  object$call$xreg <- NULL
  
  newdata %>%
    mutate(mean = biasadj(invert_transformation(object%@%"transformation"), fc$se^2)(fc$pred),
           distribution = new_fcdist(qnorm, fc$pred, sd = fc$se,
                                     transformation = invert_transformation(object%@%"transformation"),
                                     abbr = "N")
    )
}

#' @importFrom dplyr case_when
#' @export
model_sum.RW <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  drift <- length(coef(x) == 1) && "drift" == names(x$coef)
  result <- case_when(
    order[2]==1 && drift ~ "RW",
    order[2]==1 ~ "NAIVE",
    order[5]==1 ~ "SNAIVE",
    TRUE ~ "RW"
    )
  if(length(coef(x) > 1)){
    if (drift) {
      result <- paste(result, "w/ drift")
    } else {
      result <- paste("Regression with", result, "errors")
    }
  }
  result
}