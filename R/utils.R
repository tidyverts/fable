#' @import rlang
#' @importFrom purrr map
NULL

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom tsibble as_tsibble
#' @export
tsibble::as_tsibble

#' Retry with backup function
#' 
#' Attempts functions sequentially until one evaluates without error.
#' 
#' @param ... A list of functions, formulas, or atomic vectors.
#'
#'   If a __function__, it is used as is.
#'
#'   If a __formula__, e.g. `~ .x + 2`, it is converted to a function. There
#'   are three ways to refer to the arguments:
#'
#'   * For a single argument function, use `.`
#'   * For a two argument function, use `.x` and `.y`
#'   * For more arguments, use `..1`, `..2`, `..3` etc
#'
#'   This syntax allows you to create very compact anonymous functions.
#'
#'   If __character vector__, __numeric vector__, or __list__, it
#'   is converted to an extractor function. Character vectors index by name
#'   and numeric vectors index by position; use a list to index by position
#'   and name at different levels. Within a list, wrap strings in [get-attr()]
#'   to extract named attributes. If a component is not present, the value of
#'   `.default` will be returned.
#'   
#' @inheritParams purrr::possibly
#' 
#' @importFrom purrr as_mapper
#' @export
retry <- function(..., otherwise, quiet = TRUE) {
  .f <- list(...)
  
  missing_otherwise <- missing(otherwise)
  
  fn <- function(..., .list_fn_pos = 1) {
    tryCatch(as_mapper(.f[[.list_fn_pos]])(...),
             error = function(e) {
               if (!quiet)
                 message("Error: ", e$message)
               if(.list_fn_pos >= length(.f)){
                 if(missing_otherwise){
                   stop(e)
                 }
                 force(otherwise)
                 otherwise
               }
               else{
                 fn(..., .list_fn_pos = .list_fn_pos + 1)
               }
             },
             interrupt = function(e) {
               stop("Terminated by user", call. = FALSE)
             }
    )
  }
  fn
}

