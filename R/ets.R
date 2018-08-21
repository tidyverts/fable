#' @inherit forecast::ets
#' @param data A data frame
#' @param formula Model specification.
#' @param period Time series frequency for seasonal component.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% ETS(log(value) ~ season("A"))
#' 
#' @importFrom forecast ets
ETS <- function(data, formula, period = "smallest", ...){
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
    error = function(method = "Z", alpha = NULL, range = c(1e-04, 0.9999)){
      list(method = method, alpha = alpha, range = range)
    },
    trend = function(method = "Z", beta = NULL, range = c(1e-04, 0.9999),
                     damped = NULL, phi = NULL, phirange = c(0.8, 0.98)){
      list(method = method, beta = beta, range = range,
           damped = damped, phi = phi, phirange = phirange)
    },
    season = function(method = "Z", gamma = NULL, range = c(1e-04, 0.9999)){
      list(method = method, gamma = gamma, range = range)
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
  required_args %>% map(~ if(length(.x) > 1) {abort("Only one special of each type is allowed for ETS.")})
  
  args <- list(model = required_args %>% map(1) %>% map("method") %>% invoke("paste0", .),
               damped = required_args$trend[[1]]$damped,
               alpha = required_args$error[[1]]$alpha,
               beta = required_args$trend[[1]]$beta,
               gamma = required_args$season[[1]]$gamma,
               phi = required_args$trend[[1]]$phi,
               lower = required_args %>% map(1) %>% map(~ .x[["range"]][1]) %>% invoke("c", .) %>%
                 append(required_args$trend[[1]]$phirange[1]) %>%
                 as.numeric,
               upper = required_args %>% map(1) %>% map(~ .x[["range"]][2]) %>% invoke("c", .) %>%
                 append(required_args$trend[[1]]$phirange[2]) %>%
                 as.numeric
  )
  model_inputs$specials <- args
  
  # Output model
  wrap_ts_model("ets", data, model_inputs, period = period, cl = cl, ...)
}

#' @export
model_sum.ets <- function(x){
  x$method
}