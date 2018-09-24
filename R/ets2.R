#' @inherit forecast::ets
#' @param data A data frame
#' @param formula Model specification.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% ETS2(log(value) ~ season("A"))
ETS2 <- function(data, formula, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Coerce data
  data <- as_tsibble(data)
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    error = function(method = c("A", "M"), alpha = NULL, range = c(1e-04, 0.9999)){
      if (!all(is.element(method, c("A", "M")))) {
        stop("Invalid error type")
      }
      if(range[1]>range[2]){
        abort("Lower error limits must be less than upper limits")
      }
      list(method = method, alpha = alpha, range = range)
    },
    trend = function(method = c("N", "A", "Ad", "M", "Md"),
                     beta = NULL, range = c(1e-04, 0.9999),
                     phi = NULL, phirange = c(0.8, 0.98)){
      if (!all(is.element(method, c("N", "A", "Ad", "M", "Md")))) {
        stop("Invalid trend type")
      }
      if(range[1]>range[2]){
        abort("Lower trend limits must be less than upper limits")
      }
      if(phirange[1]>phirange[2]){
        abort("Lower dampening limits must be less than upper limits")
      }
      list(method = method, beta = beta, range = range,
           phi = phi, phirange = phirange)
    },
    season = function(method = c("N", "A", "M"), gamma = NULL, range = c(1e-04, 0.9999), period = "smallest"){
      if (!all(is.element(method, c("N", "A", "M")))) {
        stop("Invalid season type")
      }
      if(range[1]>range[2]){
        abort("Lower seasonal limits must be less than upper limits")
      }
      list(method = method, gamma = gamma, range = range, period = get_frequencies(period, .data))
    },
    xreg = no_xreg,
    
    .env = caller_env(),
    .required_specials = c("error", "trend", "season"),
    .vals = list(.data = data)
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  # Rebuild `ets` arguments
  parsed_args <- model_inputs$specials
  required_args <- parsed_args[c("error", "trend", "season")]
  required_args %>% map(function(.x){if(length(.x) > 1) {abort("Only one special of each type is allowed for ETS.")}})
  
  # Get response
  y <- eval_tidy(model_lhs(model_inputs$model), data = data)
}

#' @export
model_sum.ets <- function(x){
  x$method
}