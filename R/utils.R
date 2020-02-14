is.constant <- function(x) {
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  return(isTRUE(all.equal(x, y)))
}

assignSpecials <- function(x, env = caller_env()) {
  x %>%
    imap(function(.x, nm) {
      if (length(.x) > 1) warn(sprintf("Only one special for `%s` is allowed, defaulting to the first usage", nm))
      .x[[1]] %>%
        imap(function(.x, .y) assign(.y, .x, envir = env))
    })
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort(
      sprintf('The `%s` package must be installed to use this functionality. It can be installed with install.packages("%s")', pkg, pkg)
    )
  }
}

`%||%` <- function(x, y) if (is_null(x)) y else x
