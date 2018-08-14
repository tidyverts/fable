#' Multiple calls to a univariate model for each tsibble key
#' 
#' @param data A tsibble
#' @param cl A modelling call
#' 
#' @export
multi_univariate <- function(data, cl){
  multi_model(data, cl, syms(key_vars(data)))
}

#' Multiple calls to a model for mass modelling
#' 
#' @param data A tsibble
#' @param cl A modelling call
#' @param keys A set of keys to nest over
#' 
#' @export
multi_model <- function(data, cl, keys){
  nested_data <- data %>% 
    group_by(!!!keys) %>%
    nest(.key = ".data")
  
  # Re-evaluate cl in environment with split data
  out <- map(nested_data[[".data"]], function(x){
    eval_tidy(
      get_expr(cl),
      env = child_env(caller_env(), !!expr_text(get_expr(cl)$data) := x)
    )
  }) %>% 
    invoke("rbind", .)
  
  new_mable(bind_cols(nested_data[map_chr(keys, expr_text)], out))
}