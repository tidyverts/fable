#' Report information about an object
#' 
#' Displays the object in a suitable format for reporting.
#' 
#' @param object The object to report
#' @param ... Additional options for the reporting function
#' 
#' @export
report <- function(object, ...){
  UseMethod("report")
}