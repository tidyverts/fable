modelsplit <- function(object, ...){
  UseMethod("modelsplit")
}

modelsplit.dable <- function(object, ...){
  modelsplit(object[["decomp"]], ...)
}

#' @export
components <- function(object, ...){
  UseMethod("components")
}

components.dable <- function(object, ...){
  object[["decomp"]]
}