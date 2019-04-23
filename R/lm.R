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

  if(all(is.na(est[[measured_vars(est)]]))){
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  fit <- stats::lm(model_formula, .data, na.action = stats::na.exclude, ...)
  fitted <- predict(fit, .data)
  
  # Set up fit measures
  n <- length(fit$residuals)
  aic <- stats::extractAIC(fit)
  k <- aic[1] - 1
  smmry <- summary(fit)
  
  structure(
    list(
      model = fit,
      par = tibble(term = names(coef(fit))[!is.na(coef(fit))]%||%chr(),
                   !!!as_tibble(`colnames<-`(coef(smmry), c("estimate", "std.error", "statistic", "p.value")))),
      est = tibble(.fitted = fitted, .resid = !!residuals(fit)),
      fit = tibble(r.squared = smmry$r.squared, adj.r.squared = smmry$adj.r.squared, 
                   sigma = smmry$sigma, statistic = smmry$fstatistic[1]%||%NA,
                   p.value = possibly(stats::pf, NA)(smmry$fstatistic[1], 
                     smmry$fstatistic[2], smmry$fstatistic[3], lower.tail = FALSE),
                   df = smmry$df[1],
                   logLik = as.numeric(possibly(stats::logLik, NA)(fit)),
                   AIC = aic[2] + 2,
                   AICc = !!sym("AIC") + 2 * (k + 2) * (k + 3) / (n - k - 3),
                   BIC = !!sym("AIC") + (k + 2) * (log(n) - 2),
                   CV = mean((residuals(fit)/(1-stats::hatvalues(fit)))^2, na.rm = TRUE),
                   deviance = as.numeric(possibly(stats::deviance, NA)(fit)),
                   df.residual = as.numeric(possibly(stats::df.residual, NA)(fit))
      )
    ),
    class = "TSLM", origin = origin
  )
}

specials_tslm <- new_specials(
  common_xregs,
  xreg = model_xreg
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
  tslm_model <- new_model_class("TSLM", train = train_tslm,
                                specials = specials_tslm, origin = NULL)
  new_model_definition(tslm_model, !!enquo(formula), ...)
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
glance.TSLM <- function(x, ...){
  x$fit
}

#' @export
tidy.TSLM <- function(x, ...){
  x$par
}

#' @export
report.TSLM <- function(object, ...){
  cat(utils::capture.output(summary(object[["model"]]))[-1:-3], sep = "\n")
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
interpolate.TSLM <- function(object, new_data, ...){
  tm <- stats::terms(object$model)
  resp <- expr_text((tm%@%"variables")[[tm%@%"response" + 1]])
  missingVals <- is.na(new_data[[resp]])
  new_data[[resp]][missingVals] <- object$est$.fitted[missingVals]
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