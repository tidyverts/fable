#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' Default value for `NULL`.
#'
#' This infix function makes it easy to replace `NULL`s with a
#' default value. It's inspired by the way that Ruby's or operation (`||`)
#' works.
#'
#' @param x,y If `x` is NULL, will return `y`; otherwise returns
#'   `x`.
#' @export
#' @name null-default
#' @examples
#' 1 %||% 2
#' NULL %||% 2
`%||%` <- function(x, y) {
  if (is.null(x)) {
    y
  } else {
    x
  }
}

#' Infix attribute accessor
#'
#' @param x Object
#' @param name Attribute name
#' @export
#' @name get-attr
#' @examples
#' factor(1:3) %@% "levels"
#' mtcars %@% "class"
`%@%` <- function(x, name) attr(x, name, exact = TRUE)

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

