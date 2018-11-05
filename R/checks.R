check_gaps <- function(x){
  if (any(has_gaps(x)$.gaps)) {
    abort(sprintf("%s contains implicit gaps in time. You should check your data and convert implicit gaps into explicit missing values using `tsibble::fill_na()` if required.", deparse(substitute(x))))
  }
}