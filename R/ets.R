#' @inherit forecast::ets
#' @param data A data frame
#' @param formula Model specification.
#' @param restrict If TRUE (default), the models with infinite variance will not be allowed.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% ETS(log(value) ~ season("A"))
ETS <- function(data, formula, restrict = TRUE, ...){
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
    error = function(method = c("A", "M")){
      if (!all(is.element(method, c("A", "M")))) {
        stop("Invalid error type")
      }
      list(method = method)
    },
    trend = function(method = c("N", "A", "Ad"),
                     alpha = NULL, alpharange = c(1e-04, 0.9999),
                     beta = NULL, betarange = c(1e-04, 0.9999),
                     phi = NULL, phirange = c(0.8, 0.98)){
      if (!all(is.element(method, c("N", "A", "Ad", "M", "Md")))) {
        stop("Invalid trend type")
      }
      if(alpharange[1]>alpharange[2]){
        abort("Lower alpha limits must be less than upper limits")
      }
      if(betarange[1]>betarange[2]){
        abort("Lower beta limits must be less than upper limits")
      }
      if(phirange[1]>phirange[2]){
        abort("Lower phi limits must be less than upper limits")
      }
      list(method = method,
           alpha = alpha, alpharange = alpharange,
           beta = beta, betarange = betarange,
           phi = phi, phirange = phirange)
    },
    season = function(method = c("N", "A", "M"),
                      gamma = NULL, gammarange = c(1e-04, 0.9999),
                      period = "smallest"){
      if (!all(is.element(method, c("N", "A", "M")))) {
        abort("Invalid season type")
      }
      if(gammarange[1]>gammarange[2]){
        abort("Lower gamma limits must be less than upper limits")
      }
      
      m <- get_frequencies(period, .data)
      if (m < 1 || NROW(.data) <= m) {
        method <- "N"
      }
      if (m == 1) {
        if (!is.element("N", method)) {
          abort("Nonseasonal data")
        }
        method <- "N"
      }
      if (m > 24) {
        if (!is.element("N", method)) {
          abort("Frequency too high for seasonal period")
        } else if (length(method) > 1) {
          warn("I can't handle data with frequency greater than 24. Seasonality will be ignored.")
          method <- "N"
        }
      }
      list(method = method, gamma = gamma, gammarange = gammarange, period = m)
    },
    xreg = no_xreg,
    
    .env = caller_env(),
    .required_specials = c("error", "trend", "season"),
    .vals = list(.data = data)
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  # Rebuild `ets` arguments
  parsed_args <- model_inputs$specials
  ets_spec <- parsed_args[c("error", "trend", "season")]
  ets_spec %>% map(function(.x){if(length(.x) > 1) {abort("Only one special of each type is allowed for ETS.")}})
  ets_spec <- unlist(ets_spec, recursive = FALSE)
  
  # Get response
  y <- eval_tidy(model_lhs(model_inputs$model), data = data)
  
  # Build possible models
  model_opts <- expand.grid(errortype = ets_spec$error$method,
                            trendtype = ets_spec$trend$method,
                            seasontype = ets_spec$season$method,
                            stringsAsFactors = FALSE)
  model_opts$damped <- nchar(model_opts$trendtype) > 1
  model_opts$trendtype <- substr(model_opts$trendtype, 1, 1)
  
  # Remove bad models
  if(min(y) < 0){
    model_opts <- model_opts[model_opts$errortype != "M",]
  }
  if(restrict){
    restricted <- with(model_opts, 
                       (errortype == "A" & (trendtype == "M" | seasontype == "M")) | # AMM, AAM, AMA
                       (errortype == "M" & trendtype == "M" & seasontype == "A"))    # MMA
    model_opts <- model_opts[!restricted,]
  }
  
  # Find best model
  best <- NULL
  compare_ets <- function(errortype, trendtype, seasontype, damped){
    new <- possibly(quietly(etsmodel), NULL)(
      y, m = ets_spec$season$period,
      errortype = errortype, trendtype = trendtype, seasontype = seasontype, damped = damped,
      alpha = ets_spec$trend$alpha, alpharange = ets_spec$trend$alpharange,
      beta = ets_spec$trend$beta, betarange = ets_spec$trend$betarange,
      phi = ets_spec$trend$phi, phirange = ets_spec$trend$phirange,
      gamma = ets_spec$season$gamma, gammarange = ets_spec$season$gammarange,
      opt.crit = "lik", nmse = 3, bounds = "both", ...)
    
    if((new$aicc%||%Inf) < (best$aicc%||%Inf)){
      best <<- new
    }
    (new$aicc%||%Inf)
  }
  ic <- pmap_dbl(model_opts, compare_ets)
  
  best_spec <- model_opts[which.min(ic),]
  best$m <- ets_spec$season$period
  best$method <- with(best_spec,
                      paste("ETS(", errortype, ",", trendtype, ifelse(damped, "d", ""), ",", seasontype, ")", sep = ""))
  best$components <- as.character(best_spec)
  best$initstate <- best$states[1, ]
  np <- length(best$par)
  best$sigma2 <- sum(best$residuals^2, na.rm = TRUE) / (length(y) - np)
  
  # Output model
  mable(
    data,
    model = add_class(best, "ETS"),
    model_inputs
  )
}

#' @export
forecast.ETS <- function(object, newdata = NULL, ...){
  if(!is_regular(newdata)){
    abort("Forecasts must be regularly spaced")
  }
  
  errortype <- object$components[1]
  trendtype <- object$components[2]
  seasontype <- object$components[3]
  damped <- as.logical(object$components[4])
  
  fc_class <- if (errortype == "A" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
    f <- ets_fc_class1
  } else if (errortype == "M" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
    f <- ets_fc_class2
  } else if (errortype == "M" && trendtype != "M" && seasontype == "M") {
    f <- ets_fc_class3
  } else {
    abort("Forecasts from this ets method require bootstrapping which is not yet supported")
  }
  fc <- fc_class(h = NROW(newdata), last.state = object$states[NROW(object$states),],
                 trendtype, seasontype, damped, object$m, object$sigma2, object$par)
  
  construct_fc(newdata, fc$mu, sqrt(fc$var), new_fcdist(qnorm, fc$mu, sd = sqrt(fc$var), abbr = "N"))
}

#' @export
model_sum.ETS <- function(x){
  x$method
}

#' @export
print.ETS <- function(x, ...) {
  cat(paste(x$method, "\n\n"))
  ncoef <- length(x$initstate)
  
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", round(x$par["alpha"], 4), "\n"))
  if (x$components[2] != "N") {
    cat(paste("    beta  =", round(x$par["beta"], 4), "\n"))
  }
  if (x$components[3] != "N") {
    cat(paste("    gamma =", round(x$par["gamma"], 4), "\n"))
  }
  if (x$components[4] != "FALSE") {
    cat(paste("    phi   =", round(x$par["phi"], 4), "\n"))
  }
  
  cat("\n  Initial states:\n")
  cat(paste("    l =", round(x$initstate[1], 4), "\n"))
  if (x$components[2] != "N") {
    cat(paste("    b =", round(x$initstate[2], 4), "\n"))
  } else {
    x$initstate <- c(x$initstate[1], NA, x$initstate[2:ncoef])
    ncoef <- ncoef + 1
  }
  if (x$components[3] != "N") {
    cat("    s = ")
    if (ncoef <= 8) {
      cat(round(x$initstate[3:ncoef], 4))
    } else {
      cat(round(x$initstate[3:8], 4))
      cat("\n           ")
      cat(round(x$initstate[9:ncoef], 4))
    }
    cat("\n")
  }
  
  cat("\n  sigma:  ")
  cat(round(sqrt(x$sigma2), 4))
  if (!is.null(x$aic)) {
    stats <- c(x$aic, x$aicc, x$bic)
    names(stats) <- c("AIC", "AICc", "BIC")
    cat("\n\n")
    print(stats)
  }
}

#' @export
coef.ETS <- function(object, ...) {
  object$par
}