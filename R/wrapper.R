#' @importFrom tibble new_tibble
#' @importFrom dplyr grouped_df
#' @importFrom forecast msts
wrap_ts_model <- function(modelfn, data, parsed_model, period = "all", cl = "Call information lost", ...){
  period <- get_frequencies(period, data)

  # Fit model
  fit <- eval_tidy(call2(modelfn, expr(msts(!!model_lhs(parsed_model$model), !!period)), !!!parsed_model$specials, !!!dots_list(...)), data = data)

  # Backtransform
  fit$fitted <- invert_transformation(parsed_model$transformation)(fit$fitted)
  fit$x <- invert_transformation(parsed_model$transformation)(fit$x)
  
  # Fix components
  fit$call <- cl
  fit$series <- expr_text(model_lhs(parsed_model$model))
  # Output model
  mable(
    data,
    model = add_class(fit, "ts_model"),
    parsed_model
  )
}

#' @importFrom forecast forecast
#' @importFrom purrr map2
#' @importFrom stats qnorm time 
#' @importFrom utils tail
#' @importFrom dplyr pull
#' @export
forecast.ts_model <- function(object, data, newdata, ...){
  # Add smarter checking
  if(!is_regular(newdata)){
    abort("Forecasts must be regularly spaced")
  }
  
  h <- NROW(newdata)
  object <- rm_class(object, "ts_model")
  fc <- forecast(object, h=h, ...)
  
  # Assume normality
  se <- (fc$upper[,1] - fc$lower[,1])/qnorm(0.5 * (1 + fc$level[1] / 100))/2
  
  construct_fc(newdata, fc$mean, se, new_fcdist(qnorm, fc$mean, sd = se, abbr = "N"))
}

#' @export
model_sum.ts_model <- function(x){
  NextMethod()
}