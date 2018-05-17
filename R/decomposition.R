modelsplit <- function(object, ...){
  UseMethod("modelsplit")
}

modelsplit.dable <- function(object, ...){
  modelsplit(object[["decomp"]], ...)
}

#' Extract model or decomposition components
#' 
#' Extracts the decomposed components from an object, or the states from a state space model.
#' 
#' @param object A model or decomposition
#' 
#' @export
components <- function(object, ...){
  UseMethod("components")
}

components.dable <- function(object, ...){
  object[["decomp"]]
}