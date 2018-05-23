#' Construct hilo intervals
#'
#' @param lower,upper A numeric vector of values for lower and upper limits.
#' @param level Default `NULL` does not include 'level'. Otherwise values of
#' length 1 or as length of `lower`, expected between 0 and 100.
#'
#' @return A "hilo" object
#' 
#' @author Earo Wang & Mitchell O'Hara-Wild
#' 
#' @rdname hilo
#' 
#' @examples
#' new_hilo(lower = rnorm(10), upper = rnorm(10) + 5, level = 95L)
#'
#' @importFrom purrr pmap
#' @export
new_hilo <- function(lower, upper, level = NULL) {
  if (missing(lower) || missing(upper)) {
    abort("no default for `lower` or `upper`.")
  }
  if (any(upper < lower, na.rm = TRUE)) {
    abort("'upper' can't be lower than 'lower'.")
  }
  len <- length(lower)
  
  if(!is.null(level)){
    if (any(level < 0 | level > 100, na.rm = TRUE)) {
      abort("'level' can't be negative or greater than 100.")
    } else if (!(length(level) %in% c(1, len))) {
      abort(gettextf("'level' should be of length 1 or %d.", len))
    }
  }
  
  pmap(list(lower = lower, upper = upper, level = level), list) %>%
    enclass("hilo")
}

#' @rdname hilo
#' 
#' @param x Object to create hilo from
#' @param ... Additional parameters passed on to other methods
#' 
#' @export
hilo <- function(x, ...){
  UseMethod("hilo")
}

#' @export
hilo.default <- function(x, ...){
  abort(sprintf(
    "Objects of type `%s` are not supported by `hilo()`, you can create a custom `hilo` with `new_hilo()`",
    class(x)
  ))
}

#' Helpers for "hilo"
#'
#' @param x A "hilo" object.
#' @rdname helper
#'
#' @export
lower <- function(x) {
  stopifnot(is_hilo(x))
  x$lower
}

#' @rdname helper
#' @export
upper <- function(x) {
  stopifnot(is_hilo(x))
  x$upper
}

#' @rdname helper
#' @export
level <- function(x) {
  stopifnot(is_hilo(x))
  x$level
}

#' @rdname helper
#' @export
is_hilo <- function(x) {
  inherits(x, "hilo")
}

#' Validate whether values fall in the hilo
#'
#' @param x A numeric vector of values.
#' @param hilo A vector of hilo objects.
#'
#' @examples
#' rng <- new_hilo(lower = rnorm(10), upper = rnorm(10) + 5)
#' bt(0.2017, rng)
#'
#' @export
bt <- function(x, hilo) {
  stopifnot(is.numeric(x) || is_hilo(hilo))
  x >= hilo$lower & x <= hilo$upper
}

#' @importFrom purr map_dbl
#' @export
`$.hilo` <- function(x, name) {
  map_dbl(x, name)
}

#' @export
`[.hilo` <- function(x, ..., drop = TRUE) {
  enclass(NextMethod(), "hilo")
}

#' @export
c.hilo <- function(...) {
  dots_list(...) %>%
    map(`[`) %>%
    unlist(recursive = FALSE, use.names = FALSE) %>%
    enclass("hilo")
}

#' @export
print.hilo <- function(x, ..., digits = NULL) {
  cat(format(x, digits = digits), sep = "\n")
  invisible(x)
}

#' @export
format.hilo <- function(x, digits = NULL, ...) {
  format(compact_hilo(x, digits = digits), ...)
}

#' @export
is.na.hilo <- function(x) {
  # both lower and upper are NA's
  rowSums(is.na(matrix(c(x$lower, x$upper), ncol = 2))) == 2
}

#' @export
duplicated.hilo <- function(x, incomparables = FALSE, fromLast = FALSE, ...) {
  mat <- matrix(c(x$lower, x$upper, x$level), ncol = 3)
  duplicated(mat, incomparables = incomparables, fromLast = fromLast, ...)
}

#' @export
unique.hilo <- function(x, incomparables = FALSE, ...) {
  x[!duplicated(x, incomparables = incomparables, ...)]
}

#' @export
rep.hilo <- function(x, ...) {
  enclass(NextMethod(), "hilo")
}

#' @export
type_sum.hilo <- function(x) {
  "hilo"
}

#' @export
obj_sum.hilo <- function(x) {
  rep("hilo", length(x))
}

#' @export
is_vector_s3.hilo <- function(x) {
  TRUE
}

#' @export
pillar_shaft.hilo <- function(x, ...) {
  out <- compact_hilo(x)
  pillar::new_pillar_shaft_simple(out, align = "right", min_width = 10)
}

compact_hilo <- function(x, digits = NULL) {
  limit <- paste(
    format(x$lower, justify = "right", digits = digits),
    format(x$upper, justify = "right", digits = digits),
    sep = ", "
  )
  rng <- paste0("[", limit, "]")
  lvl <- level(x)
  if (is.null(lvl)) {
    return(rng)
  } else {
    paste0(rng, crayon::underline(lvl))
  }
}

#' @export
as.data.frame.hilo <- function(
  x, row.names = NULL, optional = FALSE, ...,
  nm = paste(deparse(substitute(x), width.cutoff = 500L), collapse = " ")
) {
  as.data.frame.vector(
    x, row.names = row.names, optional = optional, ...,
    nm = nm
  )
}