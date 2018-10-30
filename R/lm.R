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
#' USAccDeaths %>% as_tsibble %>% LM(log(value) ~ trend() + season())
#' 
#' library(tsibbledata)
#' olympic_running %>% 
#'   LM(Time ~ trend()) %>% 
#'   interpolate(olympic_running)
LM <- function(data, formula, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Check data
  stopifnot(is_tsibble(data))
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  origin <- min(data[[expr_text(tsibble::index(data))]])
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
  
  fit <- structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit)), estimate = coef(fit)),
      est = data %>% 
        transmute(!!model_lhs(model_formula),
                  .fitted = predict(fit, data),
                  .resid = !!model_lhs(model_formula) - .fitted)
    ),
    class = "LM", origin = origin
  )
  
  mable(
    data,
    model = fit,
    model_inputs
  )
}

#' @export
fitted.LM <- function(object, ...){
  select(object$est, ".fitted")
}

#' @export
residuals.LM <- function(object, ...){
  select(object$est, ".resid")
}

#' @export
augment.LM <- function(x, ...){
  x$est
}

#' @export
glance.LM <- function(x, ...){
  glance(x$model)
}

#' @export
tidy.LM <- function(x, ...){
  x$par
}

#' @importFrom stats predict
#' @export
forecast.LM <- function(object, new_data, ...){
  # Update bound values to special environment for re-evaluation
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  fc <- predict(object$model, new_data, se.fit = TRUE)
  
  construct_fc(new_data, fc$fit, fc$se.fit, new_fcdist(qnorm, fc$fit, sd = fc$se.fit, abbr = "N"))
}

#' @importFrom fablelite simulate
#' @export
simulate.LM <- function(object, new_data, ...){
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  pred <- predict(object$model, newdata = new_data)
  
  if(is.null(new_data[[".innov"]])){
    vars <- stats::deviance(object$model)/stats::df.residual(object$model)
    innov <- stats::rnorm(length(pred), sd = sqrt(vars))
  }
  
  transmute(new_data, .sim = pred + innov)
}

#' @export
interpolate.LM <- function(model, new_data, ...){
  resp <- response(model)
  missingVals <- is.na(new_data[[resp]])
  new_data[[resp]][missingVals] <- fitted(model)$.fitted[missingVals]
  new_data
}

#' @export
refit.LM <- function(object, new_data, reestimate = FALSE, ...){
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  fit <- stats::lm(formula(object$model$terms), data = new_data)
  fit$call <- object$model$call
  if(!reestimate){
    fit$coefficients <- object$model$coefficients
    fit$fitted.values <- predict(object$model, new_data)
    fit$residuals <- fit$model[,fit$terms%@%"response"] - fit$fitted.values
  }
  
  fit <- structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit)), estimate = coef(fit)),
      est = new_data %>% 
        transmute(!!model_lhs(formula(object$model$terms)),
                  .fitted = fit$fitted.values,
                  .resid = fit$residuals)
    ),
    class = "LM", origin = object%@%"origin"
  )
  
  mable(
    new_data,
    model = fit,
    object%@%"fable"
  )
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