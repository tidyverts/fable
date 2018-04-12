#' @importFrom tibble new_tibble
wrap_ts_model <- function(data, fn, model, period = "all", ...){
  period <- get_frequencies(period, data)
  
  # Fit model
  fit <- eval_tidy(call2(fn, expr(msts(!!model_lhs(model$model), !!period)), !!!dots_list(...)), data = data)
  
  # Backtransform
  fit$fitted <- model$backtransform(fit$fitted)
  fit$x <- model$backtransform(fit$x)
  
  # Output model
  new_tibble(list(x = list(data), model = list(structure(list(spec = model, model = fit), class = "ts_model"))), subclass = "mable")
}

#' @export
#' @importFrom forecast forecast
#' @importFrom dplyr mutate
forecast.mable <- function(object, ...){
  object %>%
    mutate(!!sym("forecast") := map(!!sym("model"), forecast, ...)) %>%
    new_tibble(subclass = "fable")
}

#' @importFrom forecast forecast
#' @export
forecast.ts_model <- function(object, ...){
  fc <- forecast(object$model, ...)
  fc$fitted <- object$spec$backtransform(fc$fitted)
  fc$upper <- object$spec$backtransform(fc$upper)
  fc$lower <- object$spec$backtransform(fc$lower)
  fc
}