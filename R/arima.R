#' @inherit forecast::Arima
#' @param data A data frame
#' @param formula Model specification.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>%
#'   as_tsibble %>%
#'   ARIMA(log(value) ~ pdq(0,1,1) + PDQ(0,1,1))
#' 
#' @importFrom forecast Arima
#' @importFrom stats model.frame
#' @importFrom tsibble n_keys
ARIMA <- function(data, formula, ...){
  # Capture call
  cl <- new_quosure(match.call())
  formula <- enquo(formula)
  
  # Coerce data
  data <- as_tsibble(data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    pdq = function(p = 0, d = 0, q = 0){
      list(order = c(p=p, d=d, q=q))
    },
    PDQ = function(P = 0, D = 0, Q = 0){
      list(seasonal = c(P=P, D=D, Q=Q))
    },
    trend = function(trend = TRUE){
      list(include.drift = trend)
    },
    xreg = specials_xreg,
    
    parent_env = get_env(cl)
  )
  
  # Parse model
  model <- data %>% 
    parse_model(formula, specials = specials)

  # Format args for call
  args <- model$args %>% 
    map(~ if(length(.x) > 1){stop("Only one special of each type is allowed in ARIMA")} else {.x[[1]]}) %>%
    set_names(NULL) %>%
    unlist(recursive = FALSE) %>%
    as.list # If no args are provided, unlist removes list structure

  fit <- eval_tidy(call2("Arima", expr(!!model_lhs(model$model)), !!!args), data = data)
  fit$fitted <- model$backtransform(fit$fitted)
  wrap_fc_model(data, fit)
}