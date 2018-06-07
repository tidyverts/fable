#' Create evaluation environment for specials
#' 
#' Allows extension packages to make use of the formula parsing of specials.
#' 
#' @param ... A named set of functions which used to parse formula inputs
#' @param parent_env The parent environment of the specials (to find other user objects)
#' @param required_specials The names of specials which must be provided (and if not, are included used with no inputs).
#' 
#' @export
new_specials_env <- function(..., parent_env = caller_env(), required_specials = NULL){
  child_env(parent_env, !!!dots_list(...)) %>%
    enclass(NULL, required_specials = required_specials)
}

tbl_xreg <- function(x){
  list(xreg = expr(tibble(!!!x)))
}

exprs_xreg <- function(x){
  x
}

no_xreg <- function(...){
  abort("Exogenous regressors are not supported for this model type.")
}