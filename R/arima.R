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
ARIMA <- function(data, formula, ...){
  # Capture call
  cl <- new_quosure(match.call())
  
  # Parse model
  model <- data %>% 
    parse_model(formula,
                specials = list(
                  pdq = function(p = 0, d = 0, q = 0){
                    list(order = c(p=p, d=d, q=q))
                  },
                  PDQ = function(P = 0, D = 0, Q = 0){
                    list(seasonal = c(P=P, D=D, Q=Q))
                  }
                ))

  # Evaluate specials
  args <- model$model %>%
    set_names(NULL) %>%
    map(
      function(special){
        if(length(special) > 1) stop("Only one of each type of special is allowed for ARIMA models.")
        eval_tidy(special[[1]], env = model$specials_env)
      }
    ) %>%
    unlist(recursive = FALSE)
  
  fit <- eval_tidy(call2("Arima", expr(!!f_lhs(formula)), !!!args), data = data)
  fit$fitted <- model$backtransform(fit$fitted)
  fit
}