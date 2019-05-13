globalVariables(c("p", "P", "q", "Q"))

#' @importFrom stats approx lm ts
train_arima <- function(.data, formula, specials, ic, stepwise = TRUE, 
                        greedy = TRUE, approximation = FALSE, 
                        order_constraint, ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by ARIMA.")
  }
  
  # Get args
  p <- d <- q <- P <- D <- Q <- period <- start.p <- start.q <- start.P <- start.Q <- NULL 
  assignSpecials(specials[c("pdq", "PDQ")])

  # Get response
  y <- x <- ts(.data[[measured_vars(.data)]], frequency = period)
  
  if(all(is.na(y))){
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  # Get xreg
  constant <- specials$xreg[[1]]$constant %||% c(TRUE, FALSE)
  specials$xreg[[1]] <- specials$xreg[[1]]$xreg
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
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
    
    D <- unname(feasts::unitroot_nsdiffs(x, differences = D, .period = period))
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
    
    d <- unname(feasts::unitroot_ndiffs(x, differences = d))
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
    offset <- with(stats::arima(y, order = c(0, d, 0), xreg = xreg, ...),
                   -2 * loglik - NROW(data) * log(sigma2))
  } else {
    method <- "CSS-ML"
  }
  
  # Find best model
  best <- NULL
  compare_arima <- function(p, d, q, P, D, Q, constant){
    if(constant){
      xreg <- cbind(xreg, constant = arima_constant(length(y), d, D, period))
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
    for (mod_spec in step_order)
    {
      est_ic <- do.call(compare_arima, model_opts[mod_spec,])
      if (isTRUE(is.finite(est_ic))) {
        break
      }
    }
  }
  
  if(is.null(best)){
    stop("Could not find an appropriate ARIMA model.")
  }
  
  # Compute ARMA roots
  ar <- best$model$phi
  ma <- best$model$theta
  arroot <- if(is_empty(ar) || !any(abs(ar) > 0)) cpl() else polyroot(c(1, -ar[seq_len(max(which(abs(ar) > 1e-8)))]))
  maroot <- if(is_empty(ma) || !any(abs(ma) > 0)) cpl() else polyroot(c(1, ma[seq_len(max(which(abs(ma) > 1e-8)))]))
  
  fit_coef <- coef(best)
  fit_se <- sqrt(diag(best$var.coef))
  if("constant" %in% names(fit_coef)){
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
                   std.error = fit_se%||%dbl()),
      est = tibble(
        .fitted = as.numeric(y - best$residuals), 
        .resid = as.numeric(best$residuals),
        .regression_resid = reg_resid
      ),
      fit = tibble(sigma = sqrt(best$sigma2),
                   logLik = best$loglik,
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
    
    constants <- map_lgl(dots, inherits, "numeric") 
    constant_given <- any(map_lgl(dots[constants], `%in%`, -1:1))
    
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    xreg <- model.frame(model_formula, data = self$data, na.action = stats::na.pass)
    
    list(
      constant = if(constant_given) as_logical(terms(xreg)%@%"intercept") else c(TRUE, FALSE),
      xreg = if(NCOL(xreg) == 0) NULL else xreg
    )
  },
  .required_specials = c("pdq", "PDQ")
)

#' Estimate an ARIMA model
#' 
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param stepwise Should stepwise be used?
#' @param greedy Should the stepwise search move to the next best option immediately?
#' @param approximation Should CSS be used during model selection?
#' @param order_constraint A logical predicate on the orders of `p`, `d`, `q`, 
#' `P`, `D` and `Q` to consider in the search.
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
                  approximation = FALSE, order_constraint = p + q + P + Q <= 5, ...){
  ic <- match.arg(ic)
  arima_model <- new_model_class("ARIMA", train = train_arima, 
                                 specials = specials_arima, origin = NULL,
                                 check = all_tsbl_checks)
  new_model_definition(arima_model, !!enquo(formula), ic = ic,
                       stepwise = stepwise, greedy = greedy, 
                       approximation = approximation,
                       order_constraint = enexpr(order_constraint), ...)
}


#' @export
fitted.ARIMA <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.ARIMA <- function(object, type = "innovation", ...){
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

#' @export
glance.ARIMA <- function(x, ...){
  x$fit
}

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
    "\nsigma^2 estimated as ", format(object$fit$sigma^2, digits = 4),
    ":  log likelihood=", format(round(object$fit$logLik, 2L)), "\n", sep = ""
  )
  
  cat("AIC=", format(round(object$fit$AIC, 2L)), sep = "")
  cat("   AICc=", format(round(object$fit$AICc, 2L)), sep = "")
  cat("   BIC=", format(round(object$fit$BIC, 2L)), "\n", sep = "")
}

#' @importFrom stats formula residuals
#' @export
forecast.ARIMA <- function(object, new_data = NULL, specials = NULL, 
                           bootstrap = FALSE, times = 5000, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced.")
  }
  if(bootstrap){
    abort("Bootstrapped forecasts for ARIMA are not yet implemented.")
  }
  
  specials$xreg[[1]] <- specials$xreg[[1]]$xreg
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
  if(object$spec$constant){
    intercept <- arima_constant(object$model$nobs + NROW(new_data),
                                object$spec$d, object$spec$D,
                                object$spec$period)
    xreg <- cbind(xreg, constant = intercept[object$model$nobs + seq_len(NROW(new_data))])
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
