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


key_vars.fable <- function(x){
  syms(setdiff(colnames(x), c("data", "model", "forecast")))
}