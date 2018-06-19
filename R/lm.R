#' Fit a linear model with time series components
#' 
#' @param data A data frame
#' @param formula Model specification.
#' @param ... Additional arguments passed to lm
#' 
#' @examples 
#' 
#' USAccDeaths %>% LM(log(value) ~ trend() + season())
LM <- function(data, formula, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Coerce data
  data <- as_tsibble(data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    trend = function(knots = NULL){
      origin <- min(data[[expr_text(index(data))]])
      trend(data, knots, origin) %>% as.matrix
    },
    season = function(period = "smallest"){
      season(data, period) %>% as_model_matrix
    },
    fourier = function(period = "smallest", K){
      origin <- min(data[[expr_text(index(data))]])
      fourier(data, period, K, origin) %>% as.matrix
    },
    #xreg is handled by lm
    parent_env = caller_env()
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula,
    specials = new_specials_env(trend = function(...){NULL},
                                season = function(...){NULL},
                                fourier = function(...){NULL},
                                xreg = function(...){NULL})
  )
  
  fit <- lm(as.formula(eval_tidy(model_inputs$model), specials), data, ...)
  fit$call <- cl
  
  mable(
    key_vals = as.list(data)[key_vars(data)],
    data = (data %>%
              grouped_df(key_vars(.)) %>%
              nest)$data,
    model = list(enclass(fit, "LM",
                         !!!model_inputs[c("model", "response", "transformation")]))
  )
}

model_sum.LM <- function(x){
  "LM"
}

as_model_matrix <- function(tbl){
  stats::model.matrix(~ ., data = tbl)[,-1]
}