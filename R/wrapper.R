#' @importFrom tibble new_tibble
wrap_ts_model <- function(modelfn, data, model, response, transformation, args, period = "all", cl = "Call information lost", ...){
  period <- get_frequencies(period, data)

  # Fit model
  fit <- eval_tidy(call2(modelfn, expr(msts(!!model_lhs(model), !!period)), !!!args, !!!dots_list(...)), data = data)
  
  # Backtransform
  fit$fitted <- invert_transformation(transformation)(fit$fitted)
  fit$x <- invert_transformation(transformation)(fit$x)
  
  # Fix components
  fit$call <- cl
  fit$series <- expr_text(response)
  
  # Output model
  new_tibble(list(data = list(data), model = list(enclass(fit, model = model, response = response, transformation = transformation, subclass = "ts_model"))), subclass = "mable")
}

#' @export
#' @importFrom forecast forecast
#' @importFrom dplyr mutate
forecast.mable <- function(object, ...){
  object %>%
    mutate(!!sym("forecast") := map2(!!sym("model"), !!sym("data"), forecast, ...)) %>%
    new_tibble(subclass = "fable")
}

#' @importFrom forecast forecast
#' @importFrom purrr map2
#' @importFrom stats qnorm time 
#' @importFrom utils tail
#' @export
forecast.ts_model <- function(object, data, bootstrap = FALSE, ...){
  if(bootstrap){
    abort("Bootstrap forecast intervals not yet supported for this model")
  }
  class(object) <- class(object)[-match("ts_model", class(object))]
  fc <- forecast(object, ...)
  # Assume normality
  se <- (fc$mean - fc$lower[,1])/qnorm(0.5 * (1 + fc$level[1] / 100))
  
  idx <- data %>% pull(!!index(.))
  future_idx <- seq(tail(idx, 1), length.out = length(fc$mean) + 1, by = time_unit(idx)) %>% tail(-1)
  
  tsibble(!!index(data) := future_idx,
          mean = fc$mean, 
          quantile = new_quantile(qnorm, fc$mean, sd = se, transformation = invert_transformation(object%@%"transformation"), abbr = "N"),
          index = !!index(data))
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