#' @importFrom purrr safely
guess_response <- function(.data){
  all_vars <- custom_error(
    measured_vars,
    "This model function does not support automatic selection of response variables. Please specify this in the model formula."
  )(.data)
  
  if(length(all_vars)!=1){
    abort("Could not automatically determine the response variable, please provide the response variable in the model specification")
  }
  
  out <- sym(all_vars[[1]])
  inform(sprintf(
    "Model not specified, defaulting to automatic modelling of the `%s` variable. Override this using the model formula.",
    expr_text(out)
  ))
  out
}