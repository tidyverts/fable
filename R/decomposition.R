modelsplit <- function(object, ...){
  UseMethod("modelsplit")
}

modelsplit.dable <- function(object, ...){
  modelsplit(object[["decomp"]], ...)
}