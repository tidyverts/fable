#' Estimate an ARIMA model
#' @param data A tsibble
#' @param formula Model specification.
#' @param ic The information criterion used to choose the best model
#' @param test The unit root test used for selecting non-seasonal differences
#' @param seasonal.test The seasonal unit root test for selecting seasonal differences
#' @param ... Further arguments for arima
#' 
#' @export
#' 
#' @examples 
#' # Manual ARIMA specification
#' USAccDeaths %>% ARIMA2(log(value) ~ pdq(0,1,1) + PDQ(0,1,1))
#' 
#' # Automatic ARIMA specification
#' tsibbledata::UKLungDeaths %>% ARIMA2(log(mdeaths) ~ pdq(0,1,1) + PDQ(0,0,1) + fdeaths + fourier(K=4))
#' 
#' @importFrom forecast Arima auto.arima
#' @importFrom stats model.matrix
#' @importFrom purrr reduce
ARIMA2 <- function(data, formula, unit_root_opts = list(), selection_opts = list(), ...){
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
  origin <- min(data[[expr_text(index(data))]])
  specials <- new_specials_env(
    !!!arima_specials,
    !!!lm_specials,
    xreg = model_xreg,
    .env = caller_env(),
    .required_specials = c("pdq", "PDQ"),
    .vals = list(.data = data, origin = origin)
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  # Get response
  y <- eval_tidy(model_lhs(model_inputs$model), data = data)
  
  # Get args
  assignSpecials(model_inputs$specials[c("pdq", "PDQ")])
  
  # Check xreg
  xreg <- model_inputs$specials[c("xreg", names(lm_specials))] %>% 
    purrr::compact() %>% 
    map(~ invoke("cbind", .x)) %>% 
    invoke("cbind", .)
  
  if(!is_empty(xreg)){
    # Check that xreg is not rank deficient
    # First check if any columns are constant
    constant_columns <- apply(xreg, 2, is.constant)
    if (any(constant_columns)) { # Remove first one
      xreg <- xreg[, -which(constant_columns)[1]]
    }
    
    # Now check if it is rank deficient
    sv <- svd(na.omit(cbind(rep(1, NROW(xreg)), xreg)))$d
    if (min(sv) / sum(sv) < .Machine$double.eps) {
      stop("xreg is rank deficient")
    }
  }
  else{
    xreg <- NULL
  }
  
  # Select differencing
  
  # Find best model
  
  # Construct appropriate output
  if(any(is.null(c(p, d, q, P, D, Q)))){
    abort("Automatic model selection is not yet supported")
  }
  
  fit <- arima(y, order = c(p, d, q),
               seasonal = list(order = c(P, D, Q), period = period),
               xreg = xreg, ...)
  
  fit$fitted <- y - fit$residuals
  
  fit$call <- cl
  
  # Output model
  mable(
    data,
    model = enclass(fit, "ARIMA2", origin = origin),
    model_inputs
  )
}

#' @export
model_sum.ARIMA2 <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  result
}

arima_specials <- list(
  pdq = function(p = NULL, d = NULL, q = NULL,
                 max.p = 5, max.d = 2, max.q = 5,
                 start.p = 2, start.q = 2){
    max.p <- min(max.p, floor(NROW(.data) / 3))
    max.q <- min(max.q, floor(NROW(.data) / 3))
    as.list(environment())
  },
  PDQ = function(P = NULL, D = NULL, Q = NULL, period = "smallest",
                 max.P = 2, max.D = 1, max.Q = 2,
                 start.P = 1, start.Q = 1){
    period <- get_frequencies(period, .data)
    max.P <- min(max.P, floor(NROW(.data) / 3 / period))
    max.Q <- min(max.Q, floor(NROW(.data) / 3 / period))
    as.list(environment())
  }
)