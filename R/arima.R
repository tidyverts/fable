globalVariables(c("p", "P", "q", "Q"))

#' @importFrom stats approx lm ts
train_arima <- function(.data, specials,  ic = "aicc",
                        stepwise = TRUE, greedy = TRUE, approximation = NULL,
                        order_constraint = p + q + P + Q <= 6, 
                        unitroot_spec = unitroot_options(), ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by ARIMA.")
  }
  
  # Get args
  p <- d <- q <- P <- D <- Q <- period <- start.p <- start.q <- start.P <- start.Q <- NULL 
  assignSpecials(specials[c("pdq", "PDQ")])

  
  # Get response
  y <- x <- ts(unclass(.data)[[measured_vars(.data)]], frequency = period)
  
  if(all(is.na(y))){
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  if(is.null(approximation)){
    approximation <- (length(x) > 150) || (period > 12)
  }
  
  # Get xreg
  constant <- specials$xreg[[1]]$constant %||% c(TRUE, FALSE)
  xreg <- specials$xreg[[1]]$xreg
  
  # Check xreg
  if(!is_empty(xreg)){
    xreg <- as.matrix(xreg)
    # Check that xreg is not rank deficient
    # First check if any columns are constant
    constant_columns <- apply(xreg, 2, is.constant)
    if(all(constant_columns)){
      xreg <- NULL
    }
    else{
      if (any(constant_columns)) {
        xreg <- xreg[, -which(constant_columns)]
      }
      
      # Now check if it is rank deficient
      sv <- svd(stats::na.omit(cbind(rep(1, NROW(xreg)), xreg)))$d
      if (min(sv) / sum(sv) < .Machine$double.eps) {
        stop("xreg is rank deficient")
      }
      
      # Finally find residuals from regression in order
      # to estimate appropriate level of differencing
      j <- !is.na(x) & !is.na(rowSums(xreg))
      x[j] <- residuals(lm(x ~ xreg))
    }
  }
  else{
    xreg <- NULL
  }
  
  diff <- function(x, differences, ...){
    if(differences == 0) return(x)
    base::diff(x, differences = differences, ...)
  }
  
  # Choose seasonal differencing
  if(length(D) > 1)
  {
    require_package("feasts")
    # Valid xregs
    
    if(!is.null(xreg)){
      keep <- map_lgl(D, function(.x){
        diff_xreg <- diff(xreg, lag = period, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      D <- D[keep]
    }
    D <- unname(feasts::unitroot_nsdiffs(stats::na.contiguous(x),
                                         alpha = unitroot_spec$nsdiffs_alpha,
                                         unitroot_fn = unitroot_spec$nsdiffs_pvalue,
                                         differences = D, .period = period))
  }
  x <- diff(x, lag = period, differences = D)
  diff_xreg <- diff(xreg, lag = period, differences = D)
  if (length(d) > 1) {
    require_package("feasts")
    
    # Valid xregs
    if(!is.null(xreg)){
      keep <- map_lgl(d, function(.x){
        diff_xreg <- diff(diff_xreg, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      d <- d[keep]
    }
    
    d <- unname(feasts::unitroot_ndiffs(stats::na.contiguous(x),
                                        alpha = unitroot_spec$ndiffs_alpha,
                                        unitroot_fn = unitroot_spec$ndiffs_pvalue,
                                        differences = d))
  }
  
  # Check number of differences selected
  if (length(D) != 1) abort("Could not find appropriate number of seasonal differences.")
  if (length(d) != 1) abort("Could not find appropriate number of non-seasonal differences.")
  if (D >= 2) {
    warn("Having more than one seasonal differences is not recommended. Please consider using only one seasonal difference.")
  } else if (D + d > 2) {
    warn("Having 3 or more differencing operations is not recommended. Please consider reducing the total number of differences.")
  }
  
  if (approximation) {
    method <- "CSS"
    offset <- with(stats::arima(y, order = c(0, d, 0), xreg = xreg),
                   -2 * loglik - NROW(data) * log(sigma2))
  } else {
    method <- "CSS-ML"
  }
  
  # Find best model
  best <- NULL
  compare_arima <- function(p, d, q, P, D, Q, constant){
    if(constant){
      intercept <- arima_constant(length(y), d, D, period)
      xreg <- if(is.null(xreg)){
        matrix(intercept, dimnames = list(NULL, "constant"))
      } else {
        cbind(xreg, intercept = intercept)
      }
    }
    
    new <- wrap_arima(
      y, order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = period),
      xreg = xreg, method = method, ..., include.mean = FALSE)
    
    if(!is.null(new)){
      nstar <- length(y) - d - D * period
      npar <- length(new$coef[new$mask]) + 1
      if (method == "CSS") {
        new$aic <- offset + nstar * log(new$sigma2) + 2 * npar
      }
      if (!is.na(new$aic)) {
        new$bic <- new$aic + npar * (log(nstar) - 2)
        new$aicc <- new$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
      }
      else {
        new$aic <- new$bic <- new$aicc <- new$ic <- Inf
      }
      # Adjust residual variance to be unbiased
      new$sigma2 <- sum(new$residuals ^ 2, na.rm = TRUE) / (nstar - npar + 1)
      
      # If automatically selecting a model
      if(NROW(model_opts) > 1){
        # Check for unit roots
        minroot <- map_dbl(list(phi = -new$model$phi,
                                theta = new$model$theta),
                           function(testvec){
                             k <- abs(testvec) > 1e-8
                             if (any(k)) {
                               last.nonzero <- max(which(k))
                               testvec <- testvec[seq_len(last.nonzero)]
                               min(abs(polyroot(c(1, testvec))))
                             }
                             else{
                               2
                             }
                           }
        )
        
        if (isTRUE(min(minroot) < 1 + 1e-2)) { # Previously 1+1e-3
          new <- NULL
        } # Don't like this model
      }
    }
    if((new[[ic]]%||%Inf) < (best[[ic]]%||%Inf)){
      best <<- new
    }
    (new[[ic]]%||%Inf)
  }
  
  model_opts <- expand.grid(p = p, d = d, q = q, P = P, D = D, Q = Q, constant = constant)
  if(NROW(model_opts) > 1){
    model_opts <- filter(model_opts, !!enexpr(order_constraint), (d + D < 2) | !constant)
    if(NROW(model_opts) == 0){
      abort("There are no ARIMA models to choose from after imposing the `order_constraint`, please consider allowing more models.")
    }
    wrap_arima <- possibly(quietly(stats::arima), NULL)
  }
  else{
    wrap_arima <- stats::arima
  }
  
  if(any((model_opts$d + model_opts$D > 1) & model_opts$constant)){
    warn("Model specification induces a quadratic or higher order polynomial trend. 
This is generally discouraged, consider removing the constant or reducing the number of differences.")
  }
  constant <- unique(model_opts$constant)
    
  if(stepwise){
    # Prepare model comparison vector
    est_ic <- rep(NA_integer_, NROW(model_opts))
    best_ic <- Inf
    
    # Initial 4 models
    initial_opts <- list(start = c(start.p, d, start.q, start.P, D, start.Q, constant[1]),
                         null = c(0, d, 0, 0, D, 0, constant[1]),
                         ar = c(max(p) > 0, d, 0, max(P) > 0, D, 0, constant[1]),
                         ma = c(0, d, max(q) > 0, 0, D, max(Q) > 0, constant[1]))
    step_order <- stats::na.omit(match(initial_opts, lapply(split(model_opts, seq_len(NROW(model_opts))), as.numeric)))
    initial <- TRUE
    
    # Stepwise search
    k <- 0
    while(NROW(model_opts[step_order,]) > 0 && k < 94){
      k <- k + 1
      
      # Evaluate model
      est_ic[step_order[1]] <- do.call(compare_arima, model_opts[step_order[1],])
      
      if(greedy && !initial){
        if(update_step <- est_ic[step_order[1]] < best_ic){
          # Update best model and score
          best_ic <- est_ic[step_order[1]]
          current <- as.numeric(model_opts[step_order[1],])
        }
      }
      else{
        if(update_step <- length(step_order) == 1){
          best_ic <- min(est_ic, na.rm = TRUE)
          current <- as.numeric(model_opts[which.min(est_ic),])
        }
      }
      
      if(update_step){
        initial <- FALSE
        # Calculate new possible steps
        dist <- apply(model_opts, 1, function(x) sum((x-current)^2))
        step_order <- order(dist, model_opts$P, model_opts$Q, model_opts$p, model_opts$q)[seq_len(sum(dist <= 2))]
        step_order <- step_order[is.na(est_ic[step_order])]
      }
      else{
        # Move to next possible step
        step_order <- step_order[-1]
      }
    }
  }
  else{
    est_ic <- pmap_dbl(model_opts, compare_arima)
  }
  
  if (approximation && !is.null(best$arma)) {
    method <- "CSS-ML"
    best <- NULL
    step_order <- order(est_ic)[seq_len(sum(!is.na(est_ic)))]
    est_ic <- rep_len(Inf, length(est_ic)) # Ignore all approximate models until re-estimated
    for (mod_spec in step_order)
    {
      est_ic[mod_spec] <- do.call(compare_arima, model_opts[mod_spec,])
      if (isTRUE(is.finite(est_ic[mod_spec]))) {
        break
      }
    }
  }
  
  if(is.null(best)){
    abort("Could not find an appropriate ARIMA model.")
  }
  
  # Compute ARMA roots
  ar <- best$model$phi
  ma <- best$model$theta
  arroot <- if(is_empty(ar) || !any(abs(ar) > 0)) cpl() else polyroot(c(1, -ar[seq_len(max(which(abs(ar) > 1e-8)))]))
  maroot <- if(is_empty(ma) || !any(abs(ma) > 0)) cpl() else polyroot(c(1, ma[seq_len(max(which(abs(ma) > 1e-8)))]))
  
  fit_coef <- coef(best)
  fit_se <- sqrt(diag(best$var.coef))
  if(is_empty(fit_se)){
    fit_se <- NULL
  }
  else if(model_opts[which.min(est_ic),"constant"] && is.null(xreg)){
    fit_coef["constant"] <- fit_coef["constant"]*(1-sum(best$model$phi))
    fit_se["constant"] <- fit_se["constant"]*(1-sum(best$model$phi))
  }
  
  # Compute regression residuals
  reg_resid <- y
  if(model_opts[which.min(est_ic),"constant"]){
    xreg <- cbind(xreg, constant = arima_constant(length(y), d, D, period))
  }
  if (!is.null(xreg)) {
    reg_resid - xreg %*% as.matrix(best$coef[(sum(best$arma[1:4]) + 1):length(best$coef)])
  }
  
  # Output model
  structure(
    list(
      par = tibble(term = names(fit_coef)%||%chr(), estimate = fit_coef%||%dbl(), 
                   std.error = fit_se%||%rep(NA, length(fit_coef))) %>%
        mutate(
          statistic = !!sym("estimate") / !!sym("std.error"),
          p.value = 2 * stats::pt(abs(!!sym("statistic")),
                                  best$nobs - length(best$coef[best$mask]),
                                  lower.tail = FALSE)
        ),
      est = tibble(
        .fitted = as.numeric(y - best$residuals), 
        .resid = as.numeric(best$residuals),
        .regression_resid = reg_resid
      ),
      fit = tibble(sigma2 = best$sigma2,
                   log_lik = best$loglik,
                   AIC = best$aic, AICc = best$aicc, BIC = best$bic,
                   ar_roots = list(arroot), ma_roots = list(maroot)),
      spec = as_tibble(model_opts[which.min(est_ic),]) %>% mutate(period = period),
      model = best,
      xreg = xreg
    ),
    class = "ARIMA"
  )
}

specials_arima <- new_specials(
  pdq = function(p = 0:5, d = 0:2, q = 0:5,
                 start.p = 2, start.q = 2){
    p <- p[p <= floor(NROW(self$data) / 3)]
    q <- q[q <= floor(NROW(self$data) / 3)]
    start.p <- p[which.min(abs(p - start.p))]
    start.q <- q[which.min(abs(q - start.q))]
    as.list(environment())
  },
  PDQ = function(P = 0:2, D = 0:1, Q = 0:2, period = NULL,
                 start.P = 1, start.Q = 1){
    period <- get_frequencies(period, self$data, .auto = "smallest")
    if(period == 1){
      # Not seasonal
      P <- 0
      D <- 0
      Q <- 0
    }
    else{
      P <- P[P <= floor(NROW(self$data) / 3 / period)]
      Q <- Q[Q <= floor(NROW(self$data) / 3 / period)]
    }
    start.P <- P[which.min(abs(P - start.P))]
    start.Q <- Q[which.min(abs(Q - start.Q))]
    as.list(environment())
  },
  common_xregs,
  xreg = function(...){
    dots <- enexprs(...)
    env <- map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
    env <- if(!is_empty(env)) get_env(env[[1]]) else base_env()
    constants <- map_lgl(dots, inherits, "numeric") 
    constant_given <- any(map_lgl(dots[constants], `%in%`, -1:1))
    
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    
    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    list(
      constant = if(constant_given) as_logical(terms(xreg)%@%"intercept") else c(TRUE, FALSE),
      xreg = if(NCOL(xreg) == 0) NULL else xreg
    )
  },
  .required_specials = c("pdq", "PDQ"),
  .xreg_specials = names(common_xregs)
)

#' Estimate an ARIMA model
#' 
#' Searches through the model space specified in the specials to identify the
#' best ARIMA model which has lowest AIC, AICc or BIC value. It is implemented
#' using [`stats::arima()`] and allows ARIMA models to be used in the fable
#' framework.
#' 
#' @aliases report.ARIMA
#' 
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param stepwise Should stepwise be used?
#' @param greedy Should the stepwise search move to the next best option immediately?
#' @param approximation Should CSS (conditional sum of squares) be used during model selection? The default (`NULL`) will use the approximation if there are more than 150 observations or if the seasonal period is greater than 12.
#' @param order_constraint A logical predicate on the orders of `p`, `d`, `q`, 
#' `P`, `D` and `Q` to consider in the search.
#' @param unitroot_spec A specification of unit root tests to use in the
#' selection of `d` and `D`. See [`unitroot_options()`] for more details.
#' @param ... Further arguments for [`stats::arima()`]
#' 
#' @section Specials:
#' 
#' \subsection{pdq}{
#' The `pdq` special is used to specify non-seasonal components of the model.
#' \preformatted{
#' pdq(p = 0:5, d = 0:2, q = 0:5,
#'     start.p = 2, start.q = 2)
#' }
#'
#' \tabular{ll}{
#'   `p`        \tab The order of the non-seasonal auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `d`        \tab The order of integration for non-seasonal differencing. If multiple values are provided, one of the values will be selected via repeated KPSS tests. \cr
#'   `q`        \tab The order of the non-seasonal moving average (MA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `start.p`  \tab If `stepwise = TRUE`, `start.p` provides the initial value for `p` for the stepwise search procedure. \cr
#'   `start.q`  \tab If `stepwise = TRUE`, `start.q` provides the initial value for `q` for the stepwise search procedure.
#' }
#' }
#' 
#' \subsection{PDQ}{
#' The `PDQ` special is used to specify seasonal components of the model.
#' \preformatted{
#' PDQ(P = 0:2, D = 0:1, Q = 0:2, period = NULL,
#'     start.P = 1, start.Q = 1)
#' }
#'
#' \tabular{ll}{
#'   `P`        \tab The order of the seasonal auto-regressive (SAR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `D`        \tab The order of integration for seasonal differencing. If multiple values are provided, one of the values will be selected via repeated heuristic tests (based on strength of seasonality from an STL decomposition). \cr
#'   `Q`        \tab The order of the seasonal moving average (SMA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `period`   \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").  \cr
#'   `start.P`  \tab If `stepwise = TRUE`, `start.P` provides the initial value for `P` for the stepwise search procedure. \cr
#'   `start.Q`  \tab If `stepwise = TRUE`, `start.Q` provides the initial value for `Q` for the stepwise search procedure.
#' }
#' }
#' 
#' \subsection{xreg}{
#' Exogenous regressors can be included in an ARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#' 
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`], where including `1` will add a constant and `0` or `-1` will remove the constant. If left out, the inclusion of a constant will be determined by minimising `ic`.
#' 
#' \preformatted{
#' xreg(...)
#' }
#' 
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)
#' }
#' }
#' 
#' @examples 
#' # Manual ARIMA specification
#' USAccDeaths %>% as_tsibble %>% 
#'   model(arima = ARIMA(log(value) ~ pdq(0,1,1) + PDQ(0,1,1)))
#' 
#' # Automatic ARIMA specification
#' library(tsibble)
#' library(dplyr)
#' tsibbledata::global_economy %>% 
#'   filter(Country == "Australia") %>%
#'   model(ARIMA(log(GDP) ~ Population))
#' 
#' @importFrom stats model.matrix
#' @export
ARIMA <- function(formula, ic = c("aicc", "aic", "bic"), stepwise = TRUE, greedy = TRUE, 
                  approximation = NULL, order_constraint = p + q + P + Q <= 6, 
                  unitroot_spec = unitroot_options(), ...){
  ic <- match.arg(ic)
  arima_model <- new_model_class("ARIMA", train = train_arima, 
                                 specials = specials_arima, origin = NULL,
                                 check = all_tsbl_checks)
  new_model_definition(arima_model, !!enquo(formula), ic = ic,
                       stepwise = stepwise, greedy = greedy, 
                       approximation = approximation,
                       order_constraint = enexpr(order_constraint),
                       unitroot_spec = unitroot_spec, ...)
}


#' Extract fitted values from a fable model
#' 
#' Extracts the fitted values.
#' 
#' @inheritParams forecast.ARIMA
#' 
#' @rdname fitted
#' @export
fitted.ARIMA <- function(object, ...){
  object$est[[".fitted"]]
}

#' Extract residuals values from a fable model
#' 
#' Extracts the residuals.
#' 
#' @inheritParams forecast.ARIMA
#' @param type The type of the residuals to extract.
#' 
#' @rdname residuals
#' @export
residuals.ARIMA <- function(object, type = c("innovation", "regression"), ...){
  type <- match.arg(type)
  if(type == "innovation"){
    object$est[[".resid"]]
  }
  else if(type == "regression"){
    object$est[[".regression_resid"]]
  }
  else{
    abort(sprintf('Residuals of `type = "%s"` are not supported by ARIMA models', type))
  }
}

#' Glance an ARIMA model
#' 
#' Construct a single row summary of the ARIMA model.
#' 
#' Contains the variance of residuals (`sigma2`), the log-likelihood (`log_lik`),
#' information criterion (`AIC`, `AICc`, `BIC`) and the characteristic roots
#' (`ar_roots` and `ma_roots`). The characteristic roots can be plotted using
#' [`feasts::gg_arma()`].
#' 
#' @inheritParams generics::glance
#' 
#' @export
glance.ARIMA <- function(x, ...){
  x$fit
}

#' Tidy a fable model
#' 
#' Returns the coefficients from the model in a `tibble` format.
#' 
#' @inheritParams generics::tidy
#' 
#' @export
tidy.ARIMA <- function(x, ...){
  x$par
}

#' @export
report.ARIMA <- function(object, ...){
  par <- rbind(tidy(object)$estimate, s.e. = tidy(object)$std.error)
  colnames(par) <- tidy(object)$term
  if (NCOL(par) > 0) {
    cat("\nCoefficients:\n")
    coef <- round(par, digits = 4)
    print.default(coef, print.gap = 2)
  }
  cat(
    "\nsigma^2 estimated as ", format(object$fit$sigma2, digits = 4),
    ":  log likelihood=", format(round(object$fit$log_lik, 2L)), "\n", sep = ""
  )
  
  cat("AIC=", format(round(object$fit$AIC, 2L)), sep = "")
  cat("   AICc=", format(round(object$fit$AICc, 2L)), sep = "")
  cat("   BIC=", format(round(object$fit$BIC, 2L)), "\n", sep = "")
}

#' Forecast a model from the fable package
#' 
#' Produces forecasts from a trained model.
#' 
#' @inheritParams fablelite::forecast
#' @param specials (passed by [`fablelite::forecast.mdl_df()`]).
#' @param bootstrap If `TRUE`, then forecast distributions are computed using simulation with resampled errors.
#' @param times The number of sample paths to use in estimating the forecast distribution when `boostrap = TRUE`.
#' 
#' @importFrom stats formula residuals
#' @rdname forecast
#' @export
forecast.ARIMA <- function(object, new_data = NULL, specials = NULL, 
                           bootstrap = FALSE, times = 5000, ...){
  if(bootstrap){
    abort("Bootstrapped forecasts for ARIMA are not yet implemented.")
  }
  xreg <- specials$xreg[[1]]$xreg
  
  if(object$spec$constant){
    intercept <- arima_constant(NROW(object$est) + NROW(new_data),
      object$spec$d, object$spec$D,
      object$spec$period)[NROW(object$est) + seq_len(NROW(new_data))]
    
    xreg <- if(is.null(xreg)){
      matrix(intercept, dimnames = list(NULL, "constant"))
    } else {
      xreg <- cbind(xreg, intercept = intercept)
    }
  }
  
  # Produce predictions
  # Remove unnecessary warning for xreg
  object$model$call$xreg <- if(!is.null(xreg)){
    expr(matrix(nrow = !!NROW(object$est), ncol = !!NCOL(xreg)))
  }
  else{
    NULL
  }
  fc <- predict(object$model, n.ahead = NROW(new_data), newxreg = xreg, ...)
  object$call$xreg <- NULL
  fc$pred <- as.numeric(fc$pred)
  fc$se <- as.numeric(fc$se)
  # Output forecasts
  construct_fc(fc$pred, fc$se, dist_normal(fc$pred, fc$se))
}


#' Refit an ARIMA model
#' 
#' Applies a fitted ARIMA model to a new dataset.
#' 
#' @inheritParams forecast.ARIMA
#' @param reestimate If `TRUE`, the coefficients for the fitted model will be re-estimated to suit the new data.
#' 
#' @importFrom stats formula residuals
#' @export
refit.ARIMA <- function(object, new_data, specials = NULL, reestimate = FALSE, ...){
  # Update data for re-evaluation
  specials$pdq[[1]] <- set_names(as.list(object$spec[c("p","d","q","p","q")]),
                                 names(specials$pdq[[1]]))
  specials$PDQ[[1]] <- set_names(as.list(object$spec[c("P","D","Q","period","P","Q")]),
                                 names(specials$PDQ[[1]]))
  if(reestimate){
    return(train_arima(new_data, specials, ...))
  }
  
  out <- train_arima(new_data, specials, fixed = object$model$coef, ...)
  out$par <- object$par
  out
}

#' Interpolate missing values from a fable model
#' 
#' Applies a model specific estimation technique to predict the values of missing values in a `tsibble`, and replace them.
#' 
#' @inheritParams forecast.ARIMA
#' 
#' @importFrom stats formula residuals
#' @export
#' @export
interpolate.ARIMA <- function(object, new_data, specials, ...){
  # Get missing values
  y <- unclass(new_data)[[measured_vars(new_data)]]
  miss_val <- which(is.na(y))
  object <- refit(object, new_data, specials, ...)$model$model
  fits <- stats::KalmanSmooth(y, object)$smooth[miss_val,,drop=FALSE] %*% as.matrix(object$Z)
  
  # Update data
  i <- (miss_val-1)%%NROW(new_data) + 1
  j <- (miss_val-1)%/%NROW(new_data) + 1
  idx_pos <- match(as_string(index(new_data)), colnames(new_data))
  j <- ifelse(j>=idx_pos, j + 1, j)
  pos <- split(i, j)
  for (i in seq_along(pos)){
    new_data[[as.numeric(names(pos)[i])]][pos[[i]]] <- fits
  }
  new_data
  new_data
}

#' @export
model_sum.ARIMA <- function(x){
  out <- sprintf("ARIMA(%i,%i,%i)%s",
                 x$spec$p, x$spec$d, x$spec$q,
                 if (any(x$spec[c("P","D","Q")] > 0)) 
                   sprintf("(%i,%i,%i)[%i]",
                           x$spec$P, x$spec$D, x$spec$Q, x$spec$period)
                 else
                   ""
  )
  if (NROW(x$par) > sum(x$spec[c("p","q","P","Q","constant")])){
    out <- sprintf("LM w/ %s errors", out)
  }
  else if (x$spec$constant){
    out <- sprintf("%s w/ %s", out, 
                   if (x$spec$d + x$spec$D == 0)
                     "mean"
                   else if (x$spec$d + x$spec$D == 1)
                     "drift"
                   else
                     "poly"
    )
  }
  out
}

arima_constant <- function(n, d, D, period){
  constant <- rep(1, n)
  if(d > 0){
    constant <- stats::diffinv(constant, differences = d, xi = rep(1, d))[seq_len(n)]
  }
  if(D > 0){
    constant <- stats::diffinv(constant, lag = period, differences = D, xi = rep(1, period*D))[seq_len(n)]
  }
  constant
}

#' Options for the unit root tests for order of integration
#' 
#' @param ndiffs_alpha,nsdiffs_alpha The level for the test specified in the `pval` functions As long as `pval < alpha`, differences will be added.
#' @param ndiffs_pvalue,nsdiffs_pvalue A function (or lambda expression) which returns the probability of the . As long as `pval < alpha`, differences will be added.
#' 
#' For the function for the seasonal p-value, the seasonal period will be provided as the `.period` argument to this function.
#' A vector of data to test is available as `.` or `.x`.
#' 
#' @seealso 
#' [feasts::unitroot_ndiffs] and [feasts::unitroot_nsdiffs]
#' 
#' @export
unitroot_options <- function(ndiffs_alpha = 0.05, nsdiffs_alpha = 0.05,
                              ndiffs_pvalue = ~feasts::unitroot_kpss(.)["kpss_pvalue"],
                              nsdiffs_pvalue = ~feasts::feat_stl(., .period)[2] < 0.64){
  as.list(environment())
}