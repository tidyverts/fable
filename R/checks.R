check_gaps <- function(x) {
  if (any(tsibble::has_gaps(x)[[".gaps"]])) {
    abort(sprintf("%s contains implicit gaps in time. You should check your data and convert implicit gaps into explicit missing values using `tsibble::fill_gaps()` if required.", deparse(substitute(x))))
  }
}

check_regular <- function(x) {
  if (!is_regular(x)) {
    abort(sprintf("%s is an irregular time series, which this model does not support. You should consider if your data can be made regular, and use `tsibble::update_tsibble(%s, regular = TRUE)` if appropriate.", deparse(substitute(x)), deparse(substitute(x))))
  }
}

check_ordered <- function(x) {
  if (!is_ordered(x)) {
    abort(sprintf(
      "%s is an unordered time series. To use this model, you first must sort the data in time order using `dplyr::arrange(%s, %s)`",
      deparse(substitute(x)), paste(c(deparse(substitute(x)), key_vars(x)), collapse = ", "), index_var(x)
    ))
  }
}

all_tsbl_checks <- function(.data) {
  check_gaps(.data)
  check_regular(.data)
  check_ordered(.data)
  if (NROW(.data) == 0) {
    abort("There is no data to model. Please provide a dataset with at least one observation.")
  }
}
