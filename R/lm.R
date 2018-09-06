#' Fit a linear model with time series components
#' 
#' @param data A data frame
#' @param formula Model specification.
#' @param ... Additional arguments passed to lm
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% LM(log(value) ~ trend() + season())
LM <- function(data, formula, ...){
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
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = data, origin = origin)
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula)
  
  if(is_formula(model_inputs$model)){
    model_formula <- set_env(model_inputs$model, new_env = specials)
  }
  else{
    model_formula <- new_formula(model_inputs$model, 1, specials)
  }
  fit <- stats::lm(model_formula, data, na.action = stats::na.exclude, ...)

  fit$call <- cl
  mable(
    data,
    model = enclass(fit, "LM", origin = origin),
    model_inputs
  )
}

#' @importFrom stats predict
#' @export
forecast.LM <- function(object, newdata, ...){
  # Update bound values to special environment for re-evaluation
  attr(object$terms, ".Environment") <- new_specials_env(
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = newdata, origin = object%@%"origin")
  )
  
  fc <- predict(object, newdata, se.fit = TRUE)
  
  construct_fc(newdata, fc$fit, fc$se.fit, new_fcdist(qnorm, fc$fit, sd = fc$se.fit, abbr = "N"))
}

#' @export
interpolate.LM <- function(model, data, ...){
  resp <- model%@%"response"
  missingVals <- is.na(data[[resp]])
  data[[resp]][missingVals] <- predict(model, newdata = data)[missingVals]
  data
}

#xreg is handled by lm
lm_specials <- list(
  trend = function(knots = NULL){
    trend(.data, knots, origin) %>% as.matrix
  },
  season = function(period = "smallest"){
    season(.data, period) %>% as_model_matrix
  },
  fourier = function(period = "smallest", K){
    fourier(.data, period, K, origin) %>% as.matrix
  }
)

#' @export
model_sum.LM <- function(x){
  "LM"
}

as_model_matrix <- function(tbl){
  stats::model.matrix(~ ., data = tbl)[,-1, drop = FALSE]
}