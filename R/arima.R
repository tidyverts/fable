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
  
  # Parse Model
  model_spec <- parse_specials(!!f_rhs(formula), specials = c("pdq", "PDQ"))
  backtransform <- eval_tidy(expr(invert_transformation(!!f_lhs(formula))), data = data)
  
  # Define specials
  pdq <- function(p = 0, d = 0, q = 0){
    list(order = eval_tidy(c(p=p, d=d, q=q), env = get_env(cl)))
  }
  PDQ <- function(P = 0, D = 0, Q = 0){
    list(seasonal = eval_tidy(c(P=P, D=D, Q=Q), env = get_env(cl)))
  }
  xreg <- function(formula){
    list(xreg = model.frame(new_formula(lhs = NULL, rhs = enexpr(formula)), data = data))
  }
  
  if(!is.null(model_spec$xreg)){
    model_spec$xreg[[1]] <- call2("xreg", !!!model_spec$xreg)
  }
  
  # Evaluate specials
  args <- model_spec %>%
    set_names(NULL) %>%
    map(
      function(special){
        if(length(special) > 1) stop("Only one of each type of special is allowed for ARIMA models.")
        eval_tidy(special[[1]], env = environment())
      }
    ) %>%
    unlist(recursive = FALSE)
  
  fit <- eval_tidy(call2("Arima", expr(!!f_lhs(formula)), !!!args), data = data)
  fit$fitted <- backtransform(fit$fitted)
  fit
}