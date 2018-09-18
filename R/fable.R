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
#' @author Rob J Hyndman, Mitchell O'Hara-Wild, Earo Wang
#'
#' Maintainer: Rob.Hyndman@monash.edu
#' @keywords package
"_PACKAGE"

globalVariables(".")

#' @rawNamespace import(rlang, except = invoke)
#' @import tsibble
#' @import fablelite
#' @importFrom tibble tibble
NULL