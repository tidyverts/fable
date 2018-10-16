#' @inherit forecast::ets
#' @param data A data frame
#' @param formula Model specification.
#' @param restrict If TRUE (default), the models with infinite variance will not be allowed.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% as_tsibble %>% ETS(log(value) ~ season("A"))
ETS <- function(data, formula, restrict = TRUE, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Check data
  stopifnot(is_tsibble(data))
  
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
  idx <- data[[expr_text(index(data))]]
  
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

  fit <- structure(
    list(
      par = tibble(term = names(best$par), estimate = best$par),
      est = data %>% 
        transmute(
          !!model_lhs(model_inputs$model),
          .fitted = best$fitted,
          .resid = best$residuals
        ),
      fit = tibble(method = with(best_spec, paste("ETS(", errortype, ",", trendtype, ifelse(damped, "d", ""), ",", seasontype, ")", sep = "")),
                   period = ets_spec$season$period,
                   sigma = sqrt(sum(best$residuals^2, na.rm = TRUE) / (length(y) - length(best$par))),
                   logLik = best$loglik, AIC = best$aic, AICc = best$aicc, BIC = best$bic,
                   MSE = best$mse, AMSE = best$amse),
      states = tsibble(
        !!!set_names(list(seq(idx[[1]], by = time_unit(idx), length.out = NROW(best$states))), expr_text(index(data))),
        !!!set_names(split(best$states, col(best$states)), colnames(best$states)),
        index = expr_text(index(data))
      ),
      spec = as_tibble(best_spec)
    ),
    class = "ETS"
  )
  
  # Output model
  mable(
    data,
    model = fit,
    model_inputs
  )
}

#' @export
forecast.ETS <- function(object, new_data = NULL, simulate = FALSE, bootstrap = FALSE, times = 5000, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced")
  }
  
  errortype <- object$spec$errortype
  trendtype <- object$spec$trendtype
  seasontype <- object$spec$seasontype
  damped <- object$spec$damped
  laststate <- as.numeric(object$states[NROW(object$states),measured_vars(object$states)])
  
  fc_class <- if (errortype == "A" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
    ets_fc_class1
  } else if (errortype == "M" && trendtype %in% c("A", "N") && seasontype %in% c("N", "A")) {
    ets_fc_class2
  } else if (errortype == "M" && trendtype != "M" && seasontype == "M") {
    ets_fc_class3
  } else {
    simulate <- TRUE
  }
  if(simulate || bootstrap){
    sim <- map(seq_len(times), function(x) simulate(object, new_data, times = times, bootstrap = bootstrap)[[".sim"]]) %>% 
      transpose %>% 
      map(as.numeric)
    pred <- .C(
      "etsforecast",
      as.double(laststate),
      as.integer(object$fit$period),
      as.integer(switch(trendtype, "N" = 0, "A" = 1, "M" = 2)),
      as.integer(switch(seasontype, "N" = 0, "A" = 1, "M" = 2)),
      as.double(ifelse(damped, object$par[["estimate"]][[object$term == "phi"]], 1)),
      as.integer(NROW(new_data)),
      as.double(numeric(NROW(new_data))),
      PACKAGE = "fable"
    )[[7]]

    construct_fc(new_data, pred, map_dbl(sim, stats::sd), sample_quantile(sim))
  }
  else{
    fc <- fc_class(h = NROW(new_data),
                   last.state = laststate,
                   trendtype, seasontype, damped, object$fit$period, object$fit$sigma^2, 
                   set_names(object$par$estimate, object$par$term))
    construct_fc(new_data, fc$mu, sqrt(fc$var), new_fcdist(qnorm, fc$mu, sd = sqrt(fc$var), abbr = "N"))
  }
}

#' @export
simulate.ETS <- function(object, new_data, bootstrap = FALSE, ...){
  if(!is_regular(new_data)){
    abort("Simulation new_data must be regularly spaced")
  }
  
  start_idx <- min(new_data[[expr_text(index(new_data))]])
  start_pos <- match(start_idx, object$states[[index(object$states)]])
  
  if(is.na(start_pos)){
    abort("The first observation index of simulation data must be within the model's training set.")
  }
  
  initstate <- as.numeric(object$states[start_pos, measured_vars(object$states)])
  
  if(is.null(new_data[[".innov"]])){
    if(bootstrap){
      new_data[[".innov"]] <- sample(stats::na.omit(residuals(object) - mean(residuals(object), na.rm = TRUE)),
                                     NROW(new_data), replace = TRUE)
    }
    else{
      new_data[[".innov"]] <- stats::rnorm(NROW(new_data), sd = object$fit$sigma)
    }
  }
  
  if (object$spec$errortype == "M") {
    new_data[[".innov"]] <- pmax(-1, new_data[[".innov"]])
  }
  
  get_par <- function(par){object$par$estimate[object$par$term==par]}
  
  result <- new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := .C(
      "etssimulate",
      as.double(initstate),
      as.integer(object$fit$period),
      as.integer(switch(object$spec$errortype, "A" = 1, "M" = 2)),
      as.integer(switch(object$spec$trendtype, "N" = 0, "A" = 1, "M" = 2)),
      as.integer(switch(object$spec$seasontype, "N" = 0, "A" = 1, "M" = 2)),
      as.double(get_par("alpha")),
      as.double(ifelse(object$spec$trendtype == "N", 0, get_par("beta"))),
      as.double(ifelse(object$spec$seasontype == "N", 0, get_par("gamma"))),
      as.double(ifelse(!object$spec$damped, 1, get_par("phi"))),
      as.integer(length(!!sym(".innov"))),
      as.double(numeric(length(!!sym(".innov")))),
      as.double(!!sym(".innov")),
      PACKAGE = "fable"
    )[[11]])
    
  if (is.na(result[[".sim"]][1])) {
    stop("Problem with multiplicative damped trend")
  }
  
  result
}

#' @export
refit.ETS <- function(object, new_data, reestimate = FALSE, reinitialise = TRUE, ...){
  est_par <- function(par){
    if(any(pos <- object$par$term==par) && !reestimate){
      object$par$estimate[pos]
    } else {
      NULL
    }
  }
  
  if(!reinitialise){
    abort("Using initial paramaters is currently not implemented.")
  }
  
  y <- new_data %>% 
    transmute(
      !!model_lhs(formula(object))
    )
  idx <- y[[expr_text(index(y))]]
  
  best <- etsmodel(
    y[[measured_vars(y)]], m = object$fit$period,
    errortype = object$spec$errortype, trendtype = object$spec$trendtype,
    seasontype = object$spec$seasontype, damped = object$spec$damped,
    alpha = est_par("alpha"), beta = est_par("beta"), phi = est_par("phi"), gamma = est_par("gamma"),
    alpharange = c(1e-04, 0.9999), betarange = c(1e-04, 0.9999),
    gammarange = c(1e-04, 0.9999), phirange = c(0.8, 0.98),
    opt.crit = "lik", nmse = 3, bounds = "both")
  
  fit <- structure(
    list(
      par = tibble(term = names(best$par), estimate = best$par),
      est = mutate(y, .fitted = best$fitted, .resid = best$residuals),
      fit = tibble(method = object$fit$method,
                   period = object$fit$period,
                   sigma = sqrt(sum(best$residuals^2, na.rm = TRUE) / (NROW(y) - length(best$par))),
                   logLik = best$loglik, AIC = best$aic, AICc = best$aicc, BIC = best$bic,
                   MSE = best$mse, AMSE = best$amse),
      states = tsibble(
        !!!set_names(list(seq(idx[[1]], by = time_unit(idx), length.out = NROW(best$states))), expr_text(index(new_data))),
        !!!set_names(split(best$states, col(best$states)), colnames(best$states)),
        index = expr_text(index(new_data))
      ),
      spec = object$spec
    ),
    class = "ETS"
  )
  
  # Output model
  mable(
    new_data,
    fit,
    object%@%"fable"
  )
}

#' @export
fitted.ETS <- function(object, ...){
  select(object$est, ".fitted")
}

#' @export
residuals.ETS <- function(object, ...){
  select(object$est, ".resid")
}

#' @export
components.ETS <- function(object, ...){
  cmp <- match(c(expr_text(index(object$states)), "l", "b", "s1"), colnames(object$states))
  out <- object$states[,stats::na.exclude(cmp)]
  `colnames<-`(out, c(expr_text(index(object$states)), "level", "slope", "season"))
}

#' @export
model_sum.ETS <- function(x){
  x$fit$method
}

#' @export
print.ETS <- function(x, ...) {
  cat(paste(x$fit$method, "\n\n"))
  ncoef <- length(measured_vars(x$states))
  
  get_par <- function(par){x$par$estimate[x$par$term==par]}
  
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", format(get_par("alpha")), "\n"))
  if (x$spec$trendtype != "N") {
    cat(paste("    beta  =", format(get_par("beta")), "\n"))
  }
  if (x$spec$seasontype != "N") {
    cat(paste("    gamma =", format(get_par("gamma")), "\n"))
  }
  if (x$spec$damped) {
    cat(paste("    phi   =", format(get_par("phi")), "\n"))
  }
  
  cat("\n  Initial states:\n")
  print.data.frame(x$states[1,measured_vars(x$states)], row.names = FALSE)
  cat("\n")
  
  cat("\n  sigma:  ")
  cat(round(x$fit$sigma, 4))
  if (!is.null(x$fit$AIC)) {
    stats <- c(AIC = x$fit$AIC, AICc = x$fit$AICc, BIC = x$fit$BIC)
    cat("\n\n")
    print(stats)
  }
}

#' @export
summary.ETS <- function(object, ...) {
  print(object)
}

#' @export
coef.ETS <- function(object, ...) {
  set_names(object$par$estimate, object$par$term)
}