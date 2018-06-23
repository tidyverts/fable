guess_response <- function(.data, ...){
  UseMethod("guess_response")
}

guess_response.tbl_ts <- function(.data, ...){
  if(length(measured_vars(.data))>1){
    abort("Could not automatically determine the response variable, please provide the response variable in the model specification")
  }
  sym(measured_vars(.data)[1])
}