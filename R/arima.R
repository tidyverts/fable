#' Estimate an ARIMA model
#' @param data A tsibble
#' @param formula Model specification.
#' @param stepwise Should stepwise be used?
#' @param greedy Should the stepwise search move to the next best option immediately?
#' @param approximation Should CSS be used during model selection?
#' @param ... Further arguments for arima
#' 
#' @export
#' 
#' @examples 
#' # Manual ARIMA specification
#' USAccDeaths %>% as_tsibble %>% ARIMA(log(value) ~ pdq(0,1,1) + PDQ(0,1,1))
#' 
#' # Automatic ARIMA specification
#' tsibbledata::UKLungDeaths %>%
#'   ARIMA(log(mdeaths) ~ pdq(0,1,1) + PDQ(0,0,1) + fdeaths + fourier(K=4))
#' 
#' @importFrom stats model.matrix
ARIMA <- function(data, formula, stepwise = TRUE, greedy = TRUE, approximation = FALSE, ...){
  # Capture user call
  cl <- call_standardise(match.call())
  
  # Check data
  stopifnot(is_tsibble(data))
  check_gaps(data)
  
  formula <- validate_model(formula, data)
  
  # Handle multivariate inputs
  if(n_keys(data) > 1){
    return(multi_univariate(data, cl))
  }
  
  # Define specials
  specials <- new_specials_env(
    !!!arima_specials,
    !!!lm_specials,
    xreg = model_xreg,
    .env = caller_env(),
    .required_specials = c("pdq", "PDQ"),
    .vals = list(.data = data, origin = min(data[[expr_text(index(data))]]))
  )
  
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  # Get response
  est <- transmute(data, !!model_lhs(model_inputs$model))
  y <- est[[measured_vars(est)]]
  
  # Get args
  p <- d <- q <- P <- D <- Q <- period <- start.p <- start.d <- start.q <- start.P <- start.D <- start.Q <- NULL 
  assignSpecials(model_inputs$specials[c("pdq", "PDQ")])
  
  # Check xreg
  xreg <- model_inputs$specials[c("xreg", names(lm_specials))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
  if(!is_empty(xreg)){
    # Check that xreg is not rank deficient
    # First check if any columns are constant
    constant_columns <- apply(xreg, 2, is.constant)
    if (any(constant_columns)) { # Remove first one
      xreg <- xreg[, -which(constant_columns)[1]]
    }
    
    # Now check if it is rank deficient
    sv <- svd(stats::na.omit(cbind(rep(1, NROW(xreg)), xreg)))$d
    if (min(sv) / sum(sv) < .Machine$double.eps) {
      stop("xreg is rank deficient")
    }
  }
  else{
    xreg <- NULL
  }
  
  # Select differencing (currently done by AIC)
  if(length(d) > 1 | length(D) > 1){
    abort("Automatic selection of differencing is currently not implemented. Please specify `d` and `D` as a single value.")
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
  compare_arima <- function(p, d, q, P, D, Q){
    new <- possibly(quietly(arima), NULL)(
      y, order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = period),
      xreg = xreg, method = method, ...)
    
    if(!is.null(new)){
      nstar <- length(y) - d - D * period
      npar <- length(new$coef) + 1
      if (approximation) {
        new$aic <- offset + nstar * log(new$sigma2) + 2 * npar
      }
      
      # Adjust residual variance to be unbiased
      new$sigma2 <- sum(new$residuals ^ 2, na.rm = TRUE) / (nstar - npar + 1)
      
      # If automatically selecting a model
      if(NROW(model_opts) > 1){
        # Check for unit roots
        minroot <- map_dbl(list(phi = new$model$phi,
                                theta = new$model$theta),
                           function(testvec){
                             k <- abs(testvec) > 1e-8
                             if (any(k)) {
                               last.nonzero <- max(which(k))
                               testvec <- testvec[1:last.nonzero]
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
    
    if((new$aic%||%Inf) < (best$aic%||%Inf)){
      best <<- new
    }
    (new$aic%||%Inf)
  }
  
  model_opts <- expand.grid(p = p, d = d, q = q, P = P, D = D, Q = Q)
  if(stepwise){
    # Prepare model comparison vector
    ic <- rep(NA_integer_, NROW(model_opts))
    best_ic <- Inf
    
    # Initial 4 models
    initial_opts <- list(start = c(start.p, start.d, start.q, start.P, start.D, start.Q),
                         null = c(0, start.d, 0, 0, start.D, 0),
                         ar = c(1, start.d, 0, 1, start.D, 0),
                         ma = c(0, start.d, 1, 0, start.D, 1))
    step_order <- stats::na.omit(match(initial_opts, lapply(split(model_opts, seq_len(NROW(model_opts))), as.numeric)))
    initial <- TRUE
    
    # Stepwise search
    k <- 0
    while(NROW(model_opts[step_order,]) > 0 && k < 94){
      k <- k + 1
      
      # Evaluate model
      ic[step_order[1]] <- do.call(compare_arima, model_opts[step_order[1],])
      
      if(greedy && !initial){
        if(update_step <- ic[step_order[1]] < best_ic){
          # Update best model and score
          best_ic <- ic[step_order[1]]
          current <- as.numeric(model_opts[step_order[1],])
        }
      }
      else{
        if(update_step <- length(step_order) == 1){
          best_ic <- min(ic, na.rm = TRUE)
          current <- as.numeric(model_opts[which.min(ic),])
        }
      }
      
      if(update_step){
        initial <- FALSE
        # Calculate new possible steps
        dist <- apply(model_opts, 1, function(x) sum((x-current)^2))
        step_order <- order(dist, model_opts$P, model_opts$Q, model_opts$p, model_opts$q)[seq_len(sum(dist <= 2))]
        step_order <- step_order[is.na(ic[step_order])]
      }
      else{
        # Move to next possible step
        step_order <- step_order[-1]
      }
    }
  }
  else{
    ic <- pmap_dbl(model_opts, compare_arima)
  }
  if (approximation && !is.null(best$arma)) {
    method <- "CSS-ML"
    best <- NULL
    step_order <- order(ic)[seq_len(sum(!is.na(ic)))]
    for (mod_spec in step_order)
    {
      ic <- do.call(compare_arima, model_opts[mod_spec,])
      if (isTRUE(is.finite(ic))) {
        break
      }
    }
  }
  
  if(is.null(best)){
    stop("Could not find an appropriate ARIMA model.")
  }
  
  # Fix model call
  best$call <- cl
  
  # Output model
  mable(
    data,
    model = 
      structure(
        list(
          par = tibble(term = names(coef(best)), estimate = coef(best)),
          est = mutate(est,
                       .fitted = as.numeric(y - best$residuals),
                       .resid = as.numeric(best$residuals)
          ),
          fit = tibble(method = model_sum(best), 
                       period = period,
                       sigma = sqrt(best$sigma2),
                       logLik = best$loglik,
                       AIC = best$aic),
          model = best
        ),
        class = "ARIMA"),
    model_inputs
  )
}

#' @export
fitted.ARIMA <- function(object, ...){
  select(object$est, ".fitted")
}

#' @export
residuals.ARIMA <- function(object, ...){
  select(object$est, ".resid")
}

#' @importFrom stats formula residuals
#' @export
forecast.ARIMA <- function(object, new_data = NULL, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced")
  }
  
  # Evaluate xreg from new_data
  specials <- new_specials_env(
    !!!arima_specials,
    !!!lm_specials,
    xreg = model_xreg,
    .env = caller_env(),
    .required_specials = c("pdq", "PDQ"),
    .vals = list(.data = new_data, origin = min(object$est[[expr_text(index(object$est))]]))
  )
  
  vals <- parse_model_rhs(model_rhs(formula(object)), new_data, specials)
  xreg <- vals$specials[c("xreg", names(lm_specials))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
  # Produce predictions
  # Remove unnecessary warning for xreg
  if(!is.null(xreg)){
    object$model$call$xreg <- expr(matrix(nrow = !!NROW(object$est), ncol = !!NCOL(xreg)))
  }
  fc <- predict(object$model, n.ahead = NROW(new_data), newxreg = xreg, ...)
  object$call$xreg <- NULL
  
  # Output forecasts
  construct_fc(new_data, fc$pred, fc$se, new_fcdist(qnorm, fc$pred, sd = fc$se, abbr = "N"))
}

#' @export
model_sum.ARIMA <- function(x){
  x$fit$method
}

#' @export
model_sum.Arima <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  result
}

arima_specials <- list(
  pdq = function(p = 0:5, d = 0:2, q = 0:5,
                 start.p = 2, start.d = 0, start.q = 2){
    p <- p[p <= floor(NROW(.data) / 3)]
    q <- q[q <= floor(NROW(.data) / 3)]
    start.p <- p[which.min(abs(p - start.p))]
    start.d <- d[which.min(abs(d - start.d))]
    start.q <- q[which.min(abs(q - start.q))]
    as.list(environment())
  },
  PDQ = function(P = 0:2, D = 0:1, Q = 0:2, period = "smallest",
                 start.P = 1, start.D = 0, start.Q = 1){
    period <- get_frequencies(period, .data)
    if(period == 1){
      # Not seasonal
      P <- 0
      D <- 0
      Q <- 0
    }
    else{
      P <- P[P <= floor(NROW(.data) / 3 / period)]
      Q <- Q[Q <= floor(NROW(.data) / 3 / period)]
    }
    start.P <- P[which.min(abs(P - start.P))]
    start.D <- D[which.min(abs(D - start.D))]
    start.Q <- Q[which.min(abs(Q - start.Q))]
    as.list(environment())
  }
)