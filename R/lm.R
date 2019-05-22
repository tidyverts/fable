lm_glance_measures <- function(fit){
  # Set up fit measures
  n <- length(fit$residuals)
  aic <- stats::extractAIC(fit)
  k <- aic[1] - 1
  smmry <- summary(fit)
  
  tibble(r.squared = smmry$r.squared, adj.r.squared = smmry$adj.r.squared, 
         sigma2 = smmry$sigma^2, statistic = smmry$fstatistic[1]%||%NA,
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
}

lm_tidy_measures <- function(fit){
  tibble(term = names(coef(fit))[!is.na(coef(fit))]%||%chr(),
         !!!as_tibble(`colnames<-`(coef(summary(fit)), c("estimate", "std.error", "statistic", "p.value"))))
}

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
  
  structure(
    list(
      model = fit,
      par = lm_tidy_measures(fit),
      est = tibble(.fitted = fitted, .resid = !!residuals(fit)),
      fit = lm_glance_measures(fit)
    ),
    class = "TSLM"
  )
}

specials_tslm <- new_specials(
  common_xregs,
  xreg = model_xreg,
  .xreg_specials = names(common_xregs)
)

#' Fit a linear model with time series components
#' 
#' @param formula Model specification.
#' @param ... Additional arguments passed to lm
#' 
#' @section Specials:
#' 
#' \subsection{xreg}{
#' Exogenous regressors can be included in an ARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#' \preformatted{
#' xreg(...)
#' }
#' 
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)
#' }
#' }
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
forecast.TSLM <- function(object, new_data, specials = NULL, bootstrap = FALSE, 
                          times = 5000, ...){
  # Update data for re-evaluation
  mdl <- environment(attr(object$model$terms, ".Environment")$fourier)$self
  mdl$data <- new_data
  
  # Intervals
  if (bootstrap){ # Compute prediction intervals using simulations
    fc <- predict(object$model, new_data, se.fit = FALSE)
    sim <- map(seq_len(times), function(x){
      imitate(object, new_data, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    dist <- dist_sim(sim)
  }  else {
    fc <- predict(object$model, new_data, se.fit = TRUE)
    se <- sqrt(fc$se.fit^2 + fc$residual.scale^2)
    fc <- fc$fit
    dist <- dist_normal(fc, se)
  }
  
  construct_fc(fc, se, dist)
}

#' @importFrom fablelite imitate
#' @export
imitate.TSLM <- function(object, new_data, bootstrap = FALSE, ...){
  # Update data for re-evaluation
  mdl <- environment(attr(object$model$terms, ".Environment")$fourier)$self
  mdl$data <- new_data
  
  res <- residuals(object$model)
  pred <- predict(object$model, newdata = new_data)
  
  if(is.null(new_data[[".innov"]])){
    if(bootstrap){
      new_data[[".innov"]] <- sample(na.omit(res) - mean(res, na.rm = TRUE),
                                     NROW(new_data), replace = TRUE)
    }
    else{
      vars <- stats::deviance(object$model)/stats::df.residual(object$model)
      new_data[[".innov"]] <- stats::rnorm(length(pred), sd = sqrt(vars))
    }
  }
  
  transmute(new_data, .sim = pred + !!sym(".innov"))
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
      par = lm_tidy_measures(fit),
      est = tibble(.fitted = fitted(fit), .resid = residuals(fit)),
      fit = lm_glance_measures(fit)
    ),
    class = "TSLM"
  )
}

#' @export
model_sum.TSLM <- function(x){
  "TSLM"
}

as_model_matrix <- function(tbl){
  stats::model.matrix(~ ., data = tbl)[,-1, drop = FALSE]
}