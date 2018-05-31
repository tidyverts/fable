guess_response <- function(.data, ...){
  useMethod("guess_response")
}

guess_response.tbl_ts <- function(.data, ...){
  sym(measured_vars(.data)[1])
}