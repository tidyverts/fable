train_tslm <- function(.data, formula, specials, ...){
  # Extract raw original data
  est <- .data
  .data <- self$data
  
  if(is_formula(formula)){
    model_formula <- set_env(formula, new_env = self$specials)
  }
  else{
    model_formula <- new_formula(formula, 1, self$specials)
  }
  fit <- stats::lm(model_formula, .data, na.action = stats::na.exclude, ...)
  fitted <- predict(fit, .data)
  
  structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit))%||%chr(),
                   !!!as_tibble(`colnames<-`(coef(summary(fit)), c("estimate", "std.error", "statistic", "p.value")))),
      est = est %>% 
        mutate(.fitted = fitted,
               .resid = !!residuals(fit))
    ),
    class = "TSLM", origin = origin
  )
}

specials_tslm <- new_specials(
  common_xregs,
  xreg = model_xreg
)

tslm_model <- R6::R6Class(NULL,
                          inherit = fablelite::model_definition,
                          public = list(
                            model = "TSLM",
                            train = train_tslm,
                            specials = specials_tslm,
                            origin = NULL
                          )
)

#' Fit a linear model with time series components
#' 
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
TSLM <- function(formula, ...){
  tslm_model$new(!!enquo(formula), ...)
}

#' @export
fitted.TSLM <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.TSLM <- function(object, ...){
  object$est[[".resid"]]
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

#' @export
report.TSLM <- function(x, ...){
  cat(capture.output(summary(x[["model"]]))[-1:-3], sep = "\n")
}

#' @importFrom stats predict
#' @export
forecast.TSLM <- function(object, new_data, specials = NULL, ...){
  # Update data for re-evaluation
  mdl <- environment(attr(object$model$terms, ".Environment")$fourier)$self
  mdl$data <- new_data
  
  fc <- predict(object$model, new_data, se.fit = TRUE)
  construct_fc(fc$fit, sqrt(fc$se.fit^2 + fc$residual.scale^2), 
               dist_normal(fc$fit, sqrt(fc$se.fit^2 + fc$residual.scale^2)))
}

#' @importFrom fablelite imitate
#' @export
imitate.TSLM <- function(object, new_data, ...){
  # Update data for re-evaluation
  mdl <- environment(attr(object$model$terms, ".Environment")$fourier)$self
  mdl$data <- new_data
  
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
refit.TSLM <- function(object, new_data, specials = NULL, reestimate = FALSE, ...){
  # Update data for re-evaluation
  mdl <- environment(attr(object$model$terms, ".Environment")$fourier)$self
  mdl$data <- new_data
  
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
        transmute(!!model_lhs(list(formula = formula(object$model$terms))),
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