train_ets <- function(.data, formula, specials, opt_crit,
                      nmse, bounds, ic, restrict = TRUE, ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by ETS.")
  }
  
  # Rebuild `ets` arguments
  ets_spec <- specials[c("error", "trend", "season")]
  ets_spec %>% map(function(.x){if(length(.x) > 1) {abort("Only one special of each type is allowed for ETS.")}})
  ets_spec <- unlist(ets_spec, recursive = FALSE)
  
  # Get response
  y <- .data[[measured_vars(.data)]]
  idx <- .data[[expr_text(index(.data))]]
  
  if(any(is.na(y))){
    abort("ETS does not support missing values.")
  }
  # Build possible models
  model_opts <- expand.grid(errortype = ets_spec$error$method,
                            trendtype = ets_spec$trend$method,
                            seasontype = ets_spec$season$method,
                            stringsAsFactors = FALSE)
  model_opts$damped <- nchar(model_opts$trendtype) > 1
  model_opts$trendtype <- substr(model_opts$trendtype, 1, 1)
  
  # Remove bad models
  if(NROW(model_opts) > 1){
    if(min(y) < 0){
      model_opts <- model_opts[model_opts$errortype != "M",]
    }
    if(restrict){
      restricted <- with(model_opts, 
                         (errortype == "A" & (trendtype == "M" | seasontype == "M")) | # AMM, AAM, AMA
                           (errortype == "M" & trendtype == "M" & seasontype == "A"))    # MMA
      model_opts <- model_opts[!restricted,]
    }
  }
  
  if(NROW(model_opts) == 0){
    abort("No valid ETS models have been allowed. Consider allowing different (more stable) models, or enabling the restricted models with `restrict = FALSE`.")
  }
  
  # Find best model
  best <- NULL
  last_error <- NULL
  compare_ets <- function(errortype, trendtype, seasontype, damped){
    new <- safely(quietly(etsmodel))(
      y, m = ets_spec$season$period,
      errortype = errortype, trendtype = trendtype, seasontype = seasontype, damped = damped,
      alpha = ets_spec$trend$alpha, alpharange = ets_spec$trend$alpha_range,
      beta = ets_spec$trend$beta, betarange = ets_spec$trend$beta_range,
      phi = ets_spec$trend$phi, phirange = ets_spec$trend$phi_range,
      gamma = ets_spec$season$gamma, gammarange = ets_spec$season$gamma_range,
      opt.crit = opt_crit, nmse = nmse, bounds = bounds, ...)
    if(!is.null(new$error)){
      last_error <<- new$error
    }
    new <- new$result
    
    if((new[[ic]]%||%Inf) < (best[[ic]]%||%Inf) || is.null(best)){
      best <<- new
    }
    new[[ic]]%||%Inf
  }
  ic <- pmap_dbl(model_opts, compare_ets)
  
  if(is.null(best)){
    abort(last_error$message)
  }
  
  best_spec <- model_opts[which.min(ic),]
  best_spec$period <- ets_spec$season$period
  
  structure(
    list(
      par = tibble(term = names(best$par), estimate = best$par),
      est = .data %>% 
        mutate(
          .fitted = best$fitted,
          .resid = best$residuals
        ),
      fit = tibble(sigma = sqrt(sum(best$residuals^2, na.rm = TRUE) / (length(y) - length(best$par))),
                   logLik = best$loglik, AIC = best$aic, AICc = best$aicc, BIC = best$bic,
                   MSE = best$mse, AMSE = best$amse, MAE = best$mae),
      states = tsibble(
        !!!set_names(list(seq(idx[[1]] - time_unit(interval(.data)),
                              by = time_unit(interval(.data)),
                              length.out = NROW(best$states))), expr_text(index(.data))),
        !!!set_names(split(best$states, col(best$states)), colnames(best$states)),
        index = !!index(.data)
      ),
      spec = as_tibble(best_spec)
    ),
    class = "ETS"
  )
}

specials_ets <- new_specials(
  error = function(method = c("A", "M")){
    if (!all(is.element(method, c("A", "M")))) {
      stop("Invalid error type")
    }
    list(method = method)
  },
  trend = function(method = c("N", "A", "Ad"),
                   alpha = NULL, alpha_range = c(1e-04, 0.9999),
                   beta = NULL, beta_range = c(1e-04, 0.9999),
                   phi = NULL, phi_range = c(0.8, 0.98)){
    if (!all(is.element(method, c("N", "A", "Ad", "M", "Md")))) {
      stop("Invalid trend type")
    }
    if(alpha_range[1]>alpha_range[2]){
      abort("Lower alpha limits must be less than upper limits")
    }
    if(beta_range[1]>beta_range[2]){
      abort("Lower beta limits must be less than upper limits")
    }
    if(phi_range[1]>phi_range[2]){
      abort("Lower phi limits must be less than upper limits")
    }
    list(method = method,
         alpha = alpha, alpha_range = alpha_range,
         beta = beta, beta_range = beta_range,
         phi = phi, phi_range = phi_range)
  },
  season = function(method = c("N", "A", "M"), period = NULL,
                    gamma = NULL, gamma_range = c(1e-04, 0.9999)){
    if (!all(is.element(method, c("N", "A", "M")))) {
      abort("Invalid season type")
    }
    if(gamma_range[1]>gamma_range[2]){
      abort("Lower gamma limits must be less than upper limits")
    }
    
    m <- get_frequencies(period, self$data, .auto = "smallest")
    if (m < 1 || NROW(self$data) <= m) {
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
    list(method = method, gamma = gamma, gamma_range = gamma_range, period = m)
  },
  xreg = no_xreg,
  .required_specials = c("error", "trend", "season")
)

#' Exponential smoothing state space model
#'
#' Returns ETS model specified by the formula.
#'
#' Based on the classification of methods as described in Hyndman et al (2008).
#'
#' The methodology is fully automatic. The model is chosen automatically if not 
#' specified. This methodology performed extremely well on the M3-competition
#' data. (See Hyndman, et al, 2002, below.)
#'
#' @param formula Model specification (see "Specials" section).
#' @param opt_crit The optimization criterion. Defaults to the log-likelihood
#' `"lik"`, but can also be set to `"mse"` (Mean Square Error), `"amse"`
#' (Average MSE over first `nmse` forecast horizons), `"sigma"` (Standard 
#' deviation of residuals), or `"mae"` (Mean Absolute Error).
#' @param nmse If `opt_crit == "amse"`, `nmse` provides the number of steps for
#' average multistep MSE (`1<=nmse<=30`).
#' @param bounds Type of parameter space to impose: `"usual"` indicates
#' all parameters must lie between specified lower and upper bounds;
#' `"admissible"` indicates parameters must lie in the admissible space;
#' `"both"` (default) takes the intersection of these regions.
#' @param ic The information criterion used in selecting the model.
#' @param restrict If TRUE (default), the models with infinite variance will not
#' be allowed.
#' @param ... Other arguments
#' 
#' @section Specials:
#'
#' \subsection{error}{
#' The `error` special is used to specify the form of the error term.
#' \preformatted{
#' error(method = c("A", "M"))
#' }
#'
#' \tabular{ll}{
#'   `method`     \tab The form of the error term: either additive ("A") or multiplicative ("M").
#' }
#' }
#' 
#' \subsection{trend}{
#' The `trend` special is used to specify the form of the trend term and associated parameters.
#' \preformatted{
#' trend(method = c("N", "A", "Ad"),
#'       alpha = NULL, alpha_range = c(1e-04, 0.9999),
#'       beta = NULL, beta_range = c(1e-04, 0.9999),
#'       phi = NULL, phi_range = c(0.8, 0.98))
#' }
#'
#' \tabular{ll}{
#'   `method`     \tab The form of the trend term: either none ("N"), additive ("A"), multiplicative ("M") or damped variants ("Ad", "Md").\cr
#'   `alpha`      \tab The value of the smoothing parameter for the level. If `alpha = 0`, the level will not change over time. Conversely, if `alpha = 1` the level will update similarly to a random walk process. \cr
#'   `alpha_range` \tab If `alpha=NULL`, `alpha_range` provides bounds for the optimised value of `alpha`.\cr
#'   `beta`       \tab The value of the smoothing parameter for the slope. If `beta = 0`, the slope will not change over time. Conversely, if `beta = 1` the level will slope will have no memory of past slopes. \cr
#'   `beta_range`  \tab If `beta=NULL`, `beta_range` provides bounds for the optimised value of `beta`.\cr
#'   `phi`        \tab The value of the dampening parameter for the slope. If `phi = 0`, the slope will be dampened immediately (no slope). Conversely, if `phi = 1` the slope will not be dampened. \cr
#'   `phi_range`   \tab If `phi=NULL`, `phi_range` provides bounds for the optimised value of `phi`.
#' }
#' }
#' 
#' \subsection{season}{
#' The `season` special is used to specify the form of the seasonal term and associated parameters.
#' \preformatted{
#' season(method = c("N", "A", "M"), period = NULL,
#'        gamma = NULL, gamma_range = c(1e-04, 0.9999))
#' }
#'
#' \tabular{ll}{
#'   `method`     \tab The form of the seasonal term: either none ("N"), additive ("A") or multiplicative ("M").\cr
#'   `period`     \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").  \cr
#'   `gamma`      \tab The value of the smoothing parameter for the seasonal pattern. If `gamma = 0`, the seasonal pattern will not change over time. Conversely, if `gamma = 1` the seasonality will have no memory of past seasonal periods. \cr
#'   `gamma_range` \tab If `gamma=NULL`, `gamma_range` provides bounds for the optimised value of `gamma`.
#' }
#' }
#'
#' 
#' @author Rob J Hyndman & Mitchell O'Hara-Wild
#' @seealso \code{\link[stats]{HoltWinters}}, \code{\link{RW}},
#' \code{\link{ARIMA}}.
#' @references Hyndman, R.J., Koehler, A.B., Snyder, R.D., and Grose, S. (2002)
#' "A state space framework for automatic forecasting using exponential
#' smoothing methods", \emph{International J. Forecasting}, \bold{18}(3),
#' 439--454.
#'
#' Hyndman, R.J., Akram, Md., and Archibald, B. (2008) "The admissible
#' parameter space for exponential smoothing models". \emph{Annals of
#' Statistical Mathematics}, \bold{60}(2), 407--426.
#'
#' Hyndman, R.J., Koehler, A.B., Ord, J.K., and Snyder, R.D. (2008)
#' \emph{Forecasting with exponential smoothing: the state space approach},
#' Springer-Verlag. \url{http://www.exponentialsmoothing.net}.
#' 
#' @export
#' 
#' @examples 
#' 
#' USAccDeaths %>% as_tsibble %>% model(ETS(log(value) ~ season("A")))
ETS <- function(formula, opt_crit = c("lik", "amse", "mse", "sigma", "mae"),
                nmse = 3, bounds = c("both", "usual", "admissible"),
                ic = c("aicc", "aic", "bic"), restrict = TRUE, ...){
  opt_crit <- match.arg(opt_crit)
  bounds <- match.arg(bounds)
  ic <- match.arg(ic)
  
  ets_model <- new_model_class("ETS", train = train_ets, specials = specials_ets,
                               check = all_tsbl_checks)
  new_model_definition(ets_model, !!enquo(formula), opt_crit = opt_crit, nmse = nmse,
                       bounds = bounds, ic = ic, restrict = restrict, ...
  )
}

#' @export
forecast.ETS <- function(object, new_data, specials = NULL, simulate = FALSE, bootstrap = FALSE, times = 5000, ...){
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
    sim <- map(seq_len(times), function(x) imitate(object, new_data, times = times, bootstrap = bootstrap)[[".sim"]]) %>% 
      transpose %>% 
      map(as.numeric)
    pred <- .C(
      "etsforecast",
      as.double(laststate),
      as.integer(object$spec$period),
      as.integer(switch(trendtype, "N" = 0, "A" = 1, "M" = 2)),
      as.integer(switch(seasontype, "N" = 0, "A" = 1, "M" = 2)),
      as.double(ifelse(damped, object$par[["estimate"]][[object$term == "phi"]], 1)),
      as.integer(NROW(new_data)),
      as.double(numeric(NROW(new_data))),
      PACKAGE = "fable"
    )[[7]]

    construct_fc(pred, map_dbl(sim, stats::sd), dist_sim(sim))
  }
  else{
    fc <- fc_class(h = NROW(new_data),
                   last.state = laststate,
                   trendtype, seasontype, damped, object$spec$period, object$fit$sigma^2, 
                   set_names(object$par$estimate, object$par$term))
    construct_fc(fc$mu, sqrt(fc$var), dist_normal(fc$mu, sqrt(fc$var)))
  }
}

#' @export
imitate.ETS <- function(object, new_data, bootstrap = FALSE, ...){
  if(!is_regular(new_data)){
    abort("Simulation new_data must be regularly spaced")
  }
  
  start_idx <- min(new_data[[expr_text(index(new_data))]])
  start_pos <- match(start_idx - time_unit(interval(new_data)), object$states[[index(object$states)]])
  
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
      as.integer(object$spec$period),
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
refit.ETS <- function(object, new_data, specials = NULL, reestimate = FALSE, reinitialise = TRUE, ...){
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
      !!parse_expr(measured_vars(object$est)[1])
    )
  idx <- y[[expr_text(index(y))]]
  
  best <- etsmodel(
    y[[measured_vars(y)]], m = object$spec$period,
    errortype = object$spec$errortype, trendtype = object$spec$trendtype,
    seasontype = object$spec$seasontype, damped = object$spec$damped,
    alpha = est_par("alpha"), beta = est_par("beta"), phi = est_par("phi"), gamma = est_par("gamma"),
    alpharange = c(1e-04, 0.9999), betarange = c(1e-04, 0.9999),
    gammarange = c(1e-04, 0.9999), phirange = c(0.8, 0.98),
    opt.crit = "lik", nmse = 3, bounds = "both")
  
  structure(
    list(
      par = tibble(term = names(best$par), estimate = best$par),
      est = tibble(.fitted = best$fitted, .resid = best$residuals),
      fit = tibble(sigma = sqrt(sum(best$residuals^2, na.rm = TRUE) / (NROW(y) - length(best$par))),
                   logLik = best$loglik, AIC = best$aic, AICc = best$aicc, BIC = best$bic,
                   MSE = best$mse, AMSE = best$amse, MAE = best$mae),
      states = tsibble(
        !!!set_names(list(seq(idx[[1]] - time_unit(interval(new_data)),
                              by = time_unit(interval(new_data)),
                              length.out = NROW(best$states))), expr_text(index(new_data))),
        !!!set_names(split(best$states, col(best$states)), colnames(best$states)),
        index = !!index(new_data)
      ),
      spec = object$spec
    ),
    class = "ETS"
  )
}

#' @export
fitted.ETS <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.ETS <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
glance.ETS <- function(x, ...){
  x$fit
}

#' @export
tidy.ETS <- function(x, ...){
  x$par
}

#' @export
components.ETS <- function(object, ...){
  spec <- object$spec
  m <- spec$period
  idx <- index(object$states)
  response <- measured_vars(object$est)[[1]]
  
  cmp <- match(c(expr_text(idx), "l", "b", "s1"), colnames(object$states))
  out <- object$states[,stats::na.exclude(cmp)]
  colnames(out) <- c(expr_text(index(object$states)), "level", "slope", "season")[!is.na(cmp)]
  if(spec$seasontype != "N"){
    seasonal_init <- tsibble(
      !!expr_text(idx) := object$states[[index(object$states)]][[1]] - rev(seq_len(m-1))*time_unit(interval(object$states)),
      season = rev(as.numeric(object$states[1,paste0("s", seq_len(m-1) + 1)])),
      index = !!idx
    )
    out <- rbind(seasonal_init, out)
    seasonalities <- list(season = list(period = m, base = NA_real_))
  }
  else{
    seasonalities <- list()
  }
  
  est_vars <- object$est %>% 
    transmute(
      !!sym(response),
      remainder = !!sym(".resid")
    )
  
  out <- left_join(out, est_vars, by = expr_text(index(object$states)))
  out <- select(out, intersect(c(expr_text(idx), response, "level", "slope", "season", "remainder"), colnames(out)))
  
  eqn <- expr(lag(!!sym("level"), 1))
  if(spec$trendtype == "A"){
    if(spec$damped){
      phi <- object$par$estimate[object$par$term=="phi"]
      eqn <- expr(!!eqn + !!phi * lag(!!sym("slope"), 1))
    }
    else{
      eqn <- expr(!!eqn + lag(!!sym("slope"), 1))
    }
  } else if(spec$trendtype == "M"){
    if(spec$damped){
      phi <- object$par$estimate[object$par$term=="phi"]
      eqn <- expr(!!eqn * lag(!!sym("slope"), 1)^!!phi)
    }
    else{
      eqn <- expr(!!eqn * lag(!!sym("slope"), 1))
    }
  }
  if(spec$seasontype == "A"){
    eqn <- expr(!!eqn + lag(!!sym("season"), !!m))
  } else if(spec$seasontype == "M"){
    eqn <- expr((!!eqn) * lag(!!sym("season"), !!m))
  }
  if(spec$errortype == "A"){
    eqn <- expr(!!eqn + !!sym("remainder"))
  } else {
    eqn <- expr((!!eqn)*(1 + !!sym("remainder")))
  }
  
  fablelite::as_dable(out, resp = !!sym(response), method = model_sum(object),
                      seasons = seasonalities, aliases = list2(!!response := eqn))
}

#' @export
model_sum.ETS <- function(x){
  with(x$spec, paste("ETS(", errortype, ",", trendtype, ifelse(damped, "d", ""), ",", seasontype, ")", sep = ""))
}

#' @export
report.ETS <- function(object, ...) {
  ncoef <- length(measured_vars(object$states))
  
  get_par <- function(par){object$par$estimate[object$par$term==par]}
  
  cat("  Smoothing parameters:\n")
  cat(paste("    alpha =", format(get_par("alpha")), "\n"))
  if (object$spec$trendtype != "N") {
    cat(paste("    beta  =", format(get_par("beta")), "\n"))
  }
  if (object$spec$seasontype != "N") {
    cat(paste("    gamma =", format(get_par("gamma")), "\n"))
  }
  if (object$spec$damped) {
    cat(paste("    phi   =", format(get_par("phi")), "\n"))
  }
  
  cat("\n  Initial states:\n")
  print.data.frame(object$states[1,measured_vars(object$states)], row.names = FALSE)
  cat("\n")
  
  cat("\n  sigma:  ")
  cat(round(object$fit$sigma, 4))
  if (!is.null(object$fit$AIC)) {
    stats <- c(AIC = object$fit$AIC, AICc = object$fit$AICc, BIC = object$fit$BIC)
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