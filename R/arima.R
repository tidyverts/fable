#' @inherit forecast::Arima
#' @param data A data frame
#' @param formula Model specification.
#' @param period Time series frequency for seasonal component.
#' @param ic The information criterion used to choose the best model
#' for auto.arima
#' @param test The unit root test used for selecting non-seasonal differences
#' @param seasonal.test The seasonal unit root test for selecting seasonal differences
#' @param ... Further arguments for Arima or auto.arima
#' 
#' @export
#' 
#' @examples 
#' # Manual ARIMA specification
#' USAccDeaths %>% ARIMA(log(value) ~ pdq(0,1,1) + PDQ(0,1,1))
#' 
#' # Automatic ARIMA specification
#' tsibbledata::UKLungDeaths %>% ARIMA(mdeaths ~ fdeaths)
#' 
#' @importFrom forecast Arima auto.arima
#' @importFrom stats model.matrix
#' @importFrom purrr reduce
ARIMA <- function(data, formula, period = "smallest", 
                  ic=c("aicc", "aic", "bic"),
                  test=c("kpss", "adf", "pp"), seasonal.test=c("seas", "ocsb", "hegy", "ch"), ...){
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
    pdq = function(p = "auto", d = "auto", q = "auto",
                   max.p = 5, max.d = 2, max.q = 5,
                   start.p = 2, start.q = 2){
      list(
        auto = list(d = d, max.p = max.p, max.q = max.q, start.p = start.p, start.q = start.q),
        manual = list(order = c(p, d, q))
      )
    },
    PDQ = function(P = "auto", D = "auto", Q = "auto",
                   max.P = 2, max.D = 1, max.Q = 2,
                   start.P = 1, start.Q = 1){
      list(
        auto = list(D = D, max.P = max.P, max.Q = max.Q, start.P = start.P, start.Q = start.Q),
        manual = list(seasonal = c(P, D, Q))
      )
    },
    trend = function(trend = TRUE){
      list(include.drift = trend)
    },
    xreg = function(...){
      model_formula <- new_formula(
        lhs = NULL,
        rhs = reduce(c(0, enexprs(...)), ~ call2("+", .x, .y))
      )
      
      list(xreg = model.matrix(model_formula, data = .data))
    },
    
    .env = caller_env(),
    .required_specials = c("pdq", "PDQ"),
    .vals = list(.data = data)
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)

  args <- model_inputs$args
  
  if(any(c(args$pdq[[1]]$manual$order, args$PDQ[[1]]$manual$seasonal) == "auto")){
    if(any(c(args$pdq[[1]]$manual$order, args$PDQ[[1]]$manual$seasonal) != "auto")){
      inform("Partial automation of parameters is not yet supported, defaulting to complete auto.arima")
    }
    args$pdq[[1]]$auto$d <- ifelse(args$pdq[[1]]$manual$order[2] == "auto", NA, args$pdq[[1]]$manual$order[2])
    args$PDQ[[1]]$auto$D <- ifelse(args$PDQ[[1]]$manual$seasonal[2] == "auto", NA, args$PDQ[[1]]$manual$seasonal[2])
    args$pdq[[1]] <- args$pdq[[1]]$auto
    args$PDQ[[1]] <- args$PDQ[[1]]$auto
    args$trend <- NULL
    args$optim[[1]] <- list(ic=ic, test=test, seasonal.test=seasonal.test)
    ts_model_fn <- "auto.arima"
  }
  else{
    args$pdq[[1]] <- args$pdq[[1]]$manual
    args$PDQ[[1]] <- args$PDQ[[1]]$manual
    ts_model_fn <- "Arima"
  }
  
  model_inputs$args <- args
  
  model_inputs <- flatten_first_args(model_inputs)
  
  # Output model
  wrap_ts_model(ts_model_fn, data, model_inputs, period = period, cl=cl, ...)
}

#' @export
model_sum.ARIMA <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  result
}