is.constant <- function(x) {
  isTRUE(all(x == x[1L]))
}

assignSpecials <- function(x, env = caller_env()) {
  imap(x, function(.x, nm) {
    if (length(.x) > 1) cli::cli_warn(sprintf("Only one special for `%s` is allowed, defaulting to the first usage", nm))
    imap(.x[[1]], function(.x, .y) assign(.y, .x, envir = env))
  })
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cli::cli_abort(
      sprintf('The `%s` package must be installed to use this functionality. It can be installed with install.packages("%s")', pkg, pkg)
    )
  }
}

`%||%` <- function(x, y) if (is_null(x)) y else x

`%0%` <- function(x, y) if (length(x) == 0) y else x

# Size of one interval step in the units of the index's arithmetic
index_step <- function(idx, interval) {
  step <- default_time_units(interval)
  if (!inherits(idx, "mixtime::mixtime")) return(step)
  granules <- if (inherits(interval, "mixtime::mt_unit")) list(interval) else interval
  if (is.list(granules) && !is.object(granules) && all(map_lgl(granules, inherits, "mixtime::mt_unit"))) {
    chronon <- mixtime::chronon_glb(idx)
    n <- tryCatch(
      sum(map_dbl(granules, function(g) mixtime::chronon_cardinality(chronon, g))),
      error = function(e) NULL
    )
    if (length(n) == 1 && is.finite(n)) step <- n
  }
  as.integer(step)
}
