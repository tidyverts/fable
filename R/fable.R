#' Forecasting Functions for Tidy Time Series
#' 
#' Methods and tools for displaying and analysing univariate time series
#' forecasts in a tidy format including exponential smoothing via state space models and
#' automatic ARIMA modelling.
#'
#' \tabular{ll}{ Package: \tab fable\cr Type: \tab Package\cr License: \tab
#' GPL3\cr LazyLoad: \tab yes\cr }
#'
#' @docType package
#' @name fable-package
#' @author Rob J Hyndman
#'
#' Maintainer: Rob.Hyndman@monash.edu
#' @keywords package
NULL # Instead of "_PACKAGE" to remove inclusion of \alias{forecast}
# "_PACKAGE"

globalVariables(".")

#' @import rlang
#' @import tsibble
#' @importFrom purrr map map2 map_lgl map_chr map_dbl
NULL

#' @importFrom dplyr mutate_if
#' @export
summary.fable <- function(object, level=c(80,95)){
  suppressWarnings(
    object %>% 
      select(!!!key_vars(object), forecast) %>%
      mutate(
        forecast = map(forecast, 
                       function(fc){
                         fc %>%
                           mutate(!!!set_names(map(level, ~ expr(hilo(!!sym("distribution"), !!.x))), paste0(level, "%"))) %>%
                           select(exclude("distribution"))
                       }
        )
      ) %>%
      enclass("lst_ts", lst_col = "forecast") %>%
      unnest(forecast) %>%
      mutate_if(is.list, enclass, "hilo")
    )
}

key_vars.fable <- function(x){
  syms(setdiff(colnames(x), c("data", "model", "forecast")))
}