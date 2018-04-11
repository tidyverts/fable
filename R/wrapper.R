#' @importFrom tibble new_tibble
wrap_fc_model <- function(data, fn, model){
  # Fit model
  fit <- eval_tidy(call2(fn, expr(!!model_lhs(model$model)), !!!flatten_first_args(model$args)), data = data)
  
  # Backtransform
  fit$fitted <- model$backtransform(fit$fitted)
  fit$x <- model$backtransform(fit$x)
  
  # Output model
  new_tibble(list(x = list(data), model = list(structure(list(spec = model, model = fit), class = "forecast_model"))), subclass = "mable")
}

#' @export
#' @importFrom forecast forecast
forecast.mable <- function(object, ...){
  object %>%
    mutate(!!sym("forecast") := map(model, forecast, ...)) %>%
    new_tibble(subclass = "tidyforecast")
}

#' @importFrom forecast forecast
#' @export
forecast.forecast_model <- function(object, ...){
  fc <- forecast(object$model, ...)
  fc$fitted <- object$spec$backtransform(fc$fitted)
  fc$upper <- object$spec$backtransform(fc$upper)
  fc$lower <- object$spec$backtransform(fc$lower)
  fc
}