#' @keywords internal
"_PACKAGE"

globalVariables(".")

#' @rawNamespace import(rlang, except = invoke)
#' @import tsibble
#' @rawNamespace import(fabletools, except = c(units_since, default_time_units))
#' @import Rcpp
#' @importFrom dplyr mutate transmute filter lag left_join select
#'
#' @useDynLib fable, .registration = TRUE
NULL
