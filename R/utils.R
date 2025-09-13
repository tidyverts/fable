is.constant <- function(x) {
  x <- as.numeric(x)
  y <- rep(x[1], length(x))
  return(isTRUE(all.equal(x, y)))
}

assignSpecials <- function(x, env = caller_env()) {
  imap(x, function(.x, nm) {
    if (length(.x) > 1) cli::cli_warn("Only one special for {.code {nm}} is allowed, defaulting to the first usage")
    imap(.x[[1]], function(.x, .y) assign(.y, .x, envir = env))
  })
}

`%||%` <- function(x, y) if (is_null(x)) y else x
