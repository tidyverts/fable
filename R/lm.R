train_tslm <- function(.data, formula, specials, ...){
  # Extract raw original data
  est <- .data
  .data <- specials[["xreg"]][[1]]
  
  origin <- .data[[expr_text(index(.data))]][[1]]
  specials <- new_specials_env(
    !!!specials_tslm,
    .env = caller_env(),
    .vals = list(.data = .data, origin = origin)
  )
  
  if(is_formula(formula)){
    model_formula <- set_env(formula, new_env = specials)
  }
  else{
    model_formula <- new_formula(formula, 1, specials)
  }
  fit <- stats::lm(model_formula, .data, na.action = stats::na.exclude, ...)
  fitted <- predict(fit, .data)
  
  structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit)), estimate = coef(fit)),
      est = est %>% 
        mutate(.fitted = fitted,
               .resid = !!residuals(fit))
    ),
    class = "TSLM", origin = origin
  )
}

specials_tslm <- list(
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
#' USAccDeaths %>% as_tsibble %>% model(TSLM(log(value) ~ trend() + season()))
#' 
#' library(tsibbledata)
#' olympic_running %>% 
#'   model(TSLM(Time ~ trend())) %>% 
#'   interpolate(olympic_running)
TSLM <- fablelite::define_model(
  train = train_tslm,
  specials = new_specials_env(
    xreg = function(...){
      .data
    },
    .env = caller_env(),
    .required_specials = "xreg"
  ) # Just keep the raw data for handling by lm()
)

#' @export
fitted.TSLM <- function(object, ...){
  select(object$est, ".fitted")
}

#' @export
residuals.TSLM <- function(object, ...){
  select(object$est, ".resid")
}

#' @export
augment.TSLM <- function(x, ...){
  x$est
}

#' @export
glance.TSLM <- function(x, ...){
  glance(x$model)
}

#' @export
tidy.TSLM <- function(x, ...){
  x$par
}

#' @importFrom stats predict
#' @export
forecast.TSLM <- function(object, new_data, ...){
  # Update bound values to special environment for re-evaluation
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!specials_tslm,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  fc <- predict(object$model, new_data, se.fit = TRUE)
  
  construct_fc(fc$fit, fc$se.fit, 
               dist_normal(fc$fit, fc$se.fit))
}

#' @importFrom fablelite simulate
#' @export
simulate.TSLM <- function(object, new_data, ...){
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!specials_tslm,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  pred <- predict(object$model, newdata = new_data)
  
  if(is.null(new_data[[".innov"]])){
    vars <- stats::deviance(object$model)/stats::df.residual(object$model)
    innov <- stats::rnorm(length(pred), sd = sqrt(vars))
  }
  
  transmute(new_data, .sim = pred + innov)
}

#' @export
interpolate.TSLM <- function(model, new_data, ...){
  resp <- measured_vars(model$est)[1]
  missingVals <- is.na(new_data[[resp]])
  new_data[[resp]][missingVals] <- model$est$.fitted[missingVals]
  new_data
}

#' @export
refit.TSLM <- function(object, new_data, reestimate = FALSE, ...){
  attr(object$model$terms, ".Environment") <- new_specials_env(
    !!!specials_tslm,
    .env = caller_env(),
    .vals = list(.data = new_data, origin = object%@%"origin")
  )
  
  fit <- stats::lm(formula(object$model$terms), data = new_data)
  fit$call <- object$model$call
  if(!reestimate){
    fit$coefficients <- object$model$coefficients
    fit$fitted.values <- predict(object$model, new_data)
    fit$residuals <- fit$model[,fit$terms%@%"response"] - fit$fitted.values
  }
  
  structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit)), estimate = coef(fit)),
      est = new_data %>% 
        transmute(!!model_lhs(formula(object$model$terms)),
                  .fitted = fit$fitted.values,
                  .resid = fit$residuals)
    ),
    class = "TSLM", origin = object%@%"origin"
  )
}

#' @export
model_sum.TSLM <- function(x){
  "TSLM"
}

as_model_matrix <- function(tbl){
  stats::model.matrix(~ ., data = tbl)[,-1, drop = FALSE]
}