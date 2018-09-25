#' @inherit forecast::ets
#' @param data A data frame
#' @param formula Model specification.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% ETS2(log(value) ~ season("A"))
ETS2 <- function(data, formula, restrict = TRUE, ...){
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
    error = function(method = c("A", "M"), alpha = NULL, range = c(1e-04, 0.9999)){
      if (!all(is.element(method, c("A", "M")))) {
        stop("Invalid error type")
      }
      if(range[1]>range[2]){
        abort("Lower error limits must be less than upper limits")
      }
      list(method = method, alpha = alpha, range = range)
    },
    trend = function(method = c("N", "A", "Ad", "M", "Md"),
                     beta = NULL, range = c(1e-04, 0.9999),
                     phi = NULL, phirange = c(0.8, 0.98)){
      if (!all(is.element(method, c("N", "A", "Ad", "M", "Md")))) {
        stop("Invalid trend type")
      }
      if(range[1]>range[2]){
        abort("Lower trend limits must be less than upper limits")
      }
      if(phirange[1]>phirange[2]){
        abort("Lower dampening limits must be less than upper limits")
      }
      list(method = method, beta = beta, range = range,
           phi = phi, phirange = phirange)
    },
    season = function(method = c("N", "A", "M"), gamma = NULL, range = c(1e-04, 0.9999), period = "smallest"){
      if (!all(is.element(method, c("N", "A", "M")))) {
        abort("Invalid season type")
      }
      if(range[1]>range[2]){
        abort("Lower seasonal limits must be less than upper limits")
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
      list(method = method, gamma = gamma, range = range, period = m)
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
      alpha = ets_spec$error$alpha, alpharange = ets_spec$error$range,
      beta = ets_spec$trend$beta, betarange = ets_spec$trend$range,
      phi = ets_spec$trend$phi, phirange = ets_spec$trend$phirange,
      gamma = ets_spec$season$gamma, gammarange = ets_spec$season$range,
      opt.crit = "lik", nmse = 3, bounds = "both", ...)
    
    if((new$aic%||%Inf) < (best$aic%||%Inf)){
      best <<- new
    }
    (new$aic%||%Inf)
  }
  
  ic <- pmap_dbl(model_opts, compare_ets)
  
  best
}

#' @export
model_sum.ets <- function(x){
  x$method
}

#' @export
print.ets <- function(x, ...) {
  cat(paste(x$method, "\n\n"))
  cat(paste("Call:\n", deparse(x$call), "\n\n"))
  ncoef <- length(x$initstate)
  if (!is.null(x$lambda)) {
    cat("  Box-Cox transformation: lambda=", round(x$lambda, 4), "\n\n")
  }
  
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