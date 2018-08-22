#' Fit a linear model with time series components
#' 
#' @param data A data frame
#' @param formula Model specification.
#' @param ... Additional arguments passed to lm
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% LM(log(value) ~ trend() + season())
LM <- function(data, formula, ...){
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
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = data, origin = min(data[[expr_text(index(data))]]))
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula)
  
  if(is_formula(model_inputs$model)){
    model_formula <- set_env(model_inputs$model, new_env = specials)
  }
  else{
    model_formula <- new_formula(model_inputs$model, 1, specials)
  }
  fit <- stats::lm(model_formula, data, na.action = na.exclude, ...)

  fit$call <- cl
  mable(
    data,
    model = add_class(fit, "LM"),
    model_inputs
  )
}

#' @importFrom stats predict
#' @export
forecast.LM <- function(object, data, newdata = NULL, h=NULL, ...){
  if(is.null(newdata)){
    if(is.null(h)){
      h <- get_frequencies("all", data) %>%
        .[.>2] %>%
        min
    }
    future_idx <- data %>% pull(!!index(.)) %>% fc_idx(h)
    newdata <- tsibble(!!!set_names(list(future_idx), expr_text(index(data))), index = !!index(data))
  }

  # Update bound values to special environment for re-evaluation
  attr(object$terms, ".Environment") <- new_specials_env(
    !!!lm_specials,
    .env = caller_env(),
    .vals = list(.data = newdata, origin = min(data[[expr_text(index(data))]]))
  )
  
  fc <- predict(object, newdata, se.fit = TRUE)
  
  newdata %>%
    mutate(mean = biasadj(invert_transformation(object%@%"transformation"), fc$se.fit^2)(fc$fit),
           distribution = new_fcdist(qnorm, fc$fit, sd = fc$se.fit,
                                     transformation = invert_transformation(object%@%"transformation"),
                                     abbr = "N")
           )
}

#' @export
interpolate.LM <- function(model, data){
  data[[model%@%"response"]] <- predict(model, newdata = data)
  data
}

#xreg is handled by lm
lm_specials <- list(
  trend = function(knots = NULL){
    trend(.data, knots, origin) %>% as.matrix
  },
  season = function(period = "smallest"){
    season(.data, period) %>% as_model_matrix
  },
  fourier = function(period = "smallest", K){
    fourier(.data, period, K, origin) %>% as.matrix
  }
)

#' @export
model_sum.LM <- function(x){
  "LM"
}

as_model_matrix <- function(tbl){
  stats::model.matrix(~ ., data = tbl)[,-1, drop = FALSE]
}