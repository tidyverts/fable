#' Estimate an ARIMA model
#' @param data A tsibble
#' @param formula Model specification.
#' @param stepwise Should stepwise be used?
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
ARIMA2 <- function(data, formula, stepwise = TRUE, ...){
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
  
  # Select differencing (currently done by AIC)
  
  # Find best model
  best <- NULL
  compare_arima <- function(p, d, q, P, D, Q){
    new <- purrr::possibly(arima, NULL)(y, order = c(p, d, q),
                                        seasonal = list(order = c(P, D, Q), period = period),
                                        xreg = xreg, ...)
    if((new$aic%||%Inf) < (best$aic%||%Inf)){
      best <<- new
    }
    (new$aic%||%Inf)
  }
  
  model_opts <- expand.grid(p = p, d = d, q = q, P = P, D = D, Q = Q)
  if(stepwise){
    ic <- rep(NA_integer_, NROW(model_opts))
    
    current <- c(start.p, start.d, start.q, start.P, start.D, start.Q)
    
    initial <- list(start = current,
                    null = c(0, start.d, 0, 0, start.D, 0),
                    ar = c(1, start.d, 0, 1, start.D, 0),
                    ma = c(0, start.d, 1, 0, start.D, 1))
    
    step_order <- na.omit(match(initial, lapply(split(model_opts, seq_len(NROW(model_opts))), as.numeric)))
    
    k <- 0
    best_ic <- Inf
    while(NROW(model_opts[step_order,]) > 0 && k < 94){
      k <- k + 1
      ic[step_order[1]] <- do.call(compare_arima, c(model_opts[step_order[1],]))
      if(ic[step_order[1]] < best_ic){
        best_ic <- ic[step_order[1]]
        current <- as.numeric(model_opts[step_order[1],])
        dist <- apply(model_opts, 1, function(x) sum((x-current)^2))
        step_order <- order(dist, model_opts$P, model_opts$Q, model_opts$p, model_opts$q)[seq_len(sum(dist <= 2))]
        step_order <- step_order[is.na(ic[step_order])]
      }
      else{
        step_order <- step_order[-1]
      }
    }
  }
  else{
    purrr::pmap(model_opts, compare_arima)
  }
  
  # Construct appropriate output
  best$fitted <- y - best$residuals
  
  best$call <- cl
  
  # Output model
  mable(
    data,
    model = enclass(best, "ARIMA2", origin = origin),
    model_inputs
  )
}

#' @export
forecast.ARIMA2 <- function(object, newdata = NULL, ...){
  if(!is_regular(newdata)){
    abort("Forecasts must be regularly spaced")
  }
  
  # Evaluate xreg from newdata
  specials <- new_specials_env(
    !!!arima_specials,
    !!!lm_specials,
    xreg = model_xreg,
    .env = caller_env(),
    .required_specials = c("pdq", "PDQ"),
    .vals = list(.data = newdata, origin = object%@%"origin")
  )
  vals <- parse_model_rhs(model_rhs(formula(object)), newdata, specials)
  xreg <- vals$specials[c("xreg", names(lm_specials))] %>% 
    purrr::compact() %>% 
    map(~ invoke("cbind", .x)) %>% 
    invoke("cbind", .)
  
  # Produce predictions
  object$call$xreg <- xreg # Bypass predict.Arima NCOL check
  fc <- predict(object, n.ahead = NROW(newdata), newxreg = xreg, ...)
  object$call$xreg <- NULL
  
  # Output forecasts
  construct_fc(newdata, fc$pred, fc$se, new_fcdist(qnorm, fc$pred, sd = fc$se, abbr = "N"))
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
  pdq = function(p = 0:5, d = 0:2, q = 0:5,
                 start.p = 2, start.d = 0, start.q = 2){
    p <- p[p <= floor(NROW(.data) / 3)]
    q <- q[q <= floor(NROW(.data) / 3)]
    start.p <- p[which.min(abs(p - start.p))]
    start.d <- d[which.min(abs(d - start.d))]
    start.q <- q[which.min(abs(q - start.q))]
    as.list(environment())
  },
  PDQ = function(P = 0:2, D = 0:1, Q = 0:2, period = "smallest",
                 start.P = 1, start.D = 0, start.Q = 1){
    period <- get_frequencies(period, .data)
    P <- P[P <= floor(NROW(.data) / 3 / period)]
    Q <- Q[Q <= floor(NROW(.data) / 3 / period)]
    start.P <- P[which.min(abs(P - start.P))]
    start.D <- D[which.min(abs(D - start.D))]
    start.Q <- Q[which.min(abs(Q - start.Q))]
    as.list(environment())
  }
)