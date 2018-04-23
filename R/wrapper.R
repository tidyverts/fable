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
    mutate(!!sym("forecast") := map2(!!sym("model"), !!sym("x"), forecast, ...)) %>%
    new_tibble(subclass = "fable")
}

#' @importFrom forecast forecast
#' @importFrom purrr map2
#' @export
forecast.ts_model <- function(object, x, bootstrap = FALSE, ...){
  if(bootstrap){
    abort("Bootstrap forecast intervals not yet supported for this model")
  }
  class(object) <- class(object)[-match("ts_model", class(object))]
  fc <- forecast(object, ...)
  # Assume normality
  se <- (fc$mean - fc$lower[,1])/qnorm(0.5 * (1 + fc$level[1] / 100))
  tsibble(!!index(x) := as.numeric(time(fc$mean)),
          mean = fc$mean, 
          quantile = map2(fc$mean, se, ~ new_quantile(qnorm, invert_transformation(object%@%"transformation"), mean = .x, sd = .y)),
          index = !!index(x))
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