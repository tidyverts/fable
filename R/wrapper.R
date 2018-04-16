#' @importFrom tibble new_tibble
wrap_ts_model <- function(data, fn, model, period = "all", ...){
  period <- get_frequencies(period, data)
  
  # Fit model
  fit <- eval_tidy(call2(fn, expr(msts(!!model_lhs(model$model), !!period)), !!!dots_list(...)), data = data)
  
  # Backtransform
  fit$fitted <- (model$transformation%@%"inverse")(fit$fitted)
  fit$x <- (model$transformation%@%"inverse")(fit$x)
  
  # Fix components
  fit$series <- expr_text(model_lhs(model$model))
  
  # Output model
  new_tibble(list(x = list(data), model = list(enclass(fit, !!!model, subclass = "ts_model"))), subclass = "mable")
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
  class(object) <- class(object)[-match("ts_model", class(object))]
  fc <- forecast(object, ...)
  fc$fitted <- (object%@%"backtransform")(fc$fitted)
  fc$upper <- (object%@%"backtransform")(fc$upper)
  fc$lower <- (object%@%"backtransform")(fc$lower)
  fc
}

#' @export
type_sum.ts_model <- function(x){
  model_sum(x)
}

model_sum <- function(x){
  UseMethod("model_sum")
}

model_sum.ts_model <- function(x){
  x %>%
    rm_class("ts_model") %>%
    model_sum
}