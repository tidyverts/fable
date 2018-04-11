#' @importFrom tibble new_tibble
wrap_fc_model <- function(data, fit){
  new_tibble(list(x = list(data), model = list(fit)), subclass = "forecast_model")
}

#' @export
#' @importFrom forecast forecast
forecast.forecast_model <- function(model, ...){
  model %>%
    mutate(!!sym("forecast") := map(model, forecast, ...)) %>%
    new_tibble(subclass = "tidyforecast")
}