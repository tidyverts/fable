#' @keywords internal
"_PACKAGE"

globalVariables(".")

#' @rawNamespace import(rlang, except = invoke)
#' @import tsibble
#' @import fablelite
#' @import Rcpp
#' @importFrom dplyr mutate transmute filter lag left_join select
#' 
#' @useDynLib fable, .registration = TRUE
NULL