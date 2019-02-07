#' @importFrom stats approx lm ts
train_arima <- function(.data, formula, specials, stepwise = TRUE, 
                        greedy = TRUE, approximation = FALSE, ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by ARIMA.")
  }
  
  # Get args
  p <- d <- q <- P <- D <- Q <- period <- start.p <- start.q <- start.P <- start.Q <- NULL 
  assignSpecials(specials[c("pdq", "PDQ")])

  # Get response
  y <- x <- ts(.data[[measured_vars(.data)]], frequency = period)
  
  # Check xreg
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
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
    # Valid xregs
    if(!is.null(xreg)){
      keep <- map_lgl(D, function(.x){
        diff_xreg <- diff(xreg, lag = period, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      D <- D[keep]
    }
    # Non-missing y
    keep <- map_lgl(D, function(.x){
      dx <- diff(x, lag = period, differences = .x)
      !all(is.na(dx))
    })
    D <- D[keep]
    
    # Estimate the test
    keep <- map_lgl(D[-1]-1, function(.x) {
      seas_heuristic(diff(x, lag = period, differences = .x)) > 0.64
    })
    D <- max(D[c(TRUE, keep)], na.rm = TRUE)
  }
  x <- diff(x, lag = period, differences = D)
  diff_xreg <- diff(xreg, lag = period, differences = D)
  
  if (length(d) > 1) {
    # Valid xregs
    if(!is.null(xreg)){
      keep <- map_lgl(d, function(.x){
        diff_xreg <- diff(diff_xreg, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      d <- d[keep]
    }
    # Non-missing y
    keep <- map_lgl(d, function(.x){
      dx <- diff(x, differences = .x)
      !all(is.na(dx))
    })
    d <- d[keep]
    
    # Estimate the test
    keep <- map_lgl(d[-1]-1, function(.x) {
      test <- kpss_test(diff(x, differences = .x))
      approx(test$cval[1,], as.numeric(sub("pct", "", colnames(test$cval)))/100, xout=test$teststat[1], rule=2)$y < 0.05
    })
    d <- max(d[c(TRUE, keep)], na.rm = TRUE)
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
  compare_arima <- function(p, d, q, P, D, Q){
    new <- possibly(quietly(stats::arima), NULL)(
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
    initial_opts <- list(start = c(start.p, d, start.q, start.P, D, start.Q),
                         null = c(0, d, 0, 0, D, 0),
                         ar = c(1, d, 0, 1, D, 0),
                         ma = c(0, d, 1, 0, D, 1))
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
  
  # Output model
  structure(
    list(
      par = tibble(term = names(coef(best))%||%chr(), estimate = coef(best)%||%dbl()),
      est = mutate(.data,
                   .fitted = as.numeric(y - best$residuals),
                   .resid = as.numeric(best$residuals)
      ),
      fit = tibble(sigma = sqrt(best$sigma2),
                   logLik = best$loglik,
                   AIC = best$aic),
      spec = tibble(period = period),
      model = best
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
  PDQ = function(P = 0:2, D = 0:1, Q = 0:2, period = "smallest",
                 start.P = 1, start.Q = 1){
    period <- get_frequencies(period, self$data)
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
  xreg = model_xreg,
  .required_specials = c("pdq", "PDQ")
)

arima_model <- R6::R6Class(NULL,
                           inherit = fablelite::model_definition,
                           public = list(
                             model = "ARIMA",
                             train = train_arima,
                             specials = specials_arima,
                             origin = NULL,
                             check = function(.data){
                               check_gaps(.data)
                             }
                           )
)

#' Estimate an ARIMA model
#' 
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
#' USAccDeaths %>% as_tsibble %>% 
#'   model(arima = ARIMA(log(value) ~ pdq(0,1,1) + PDQ(0,1,1)))
#' 
#' # Automatic ARIMA specification
#' tsibbledata::UKLungDeaths %>%
#'   model(arima = ARIMA(log(mdeaths) ~ pdq(0,1,1) + PDQ(0,0,1) + 
#'                       fdeaths + fourier(K=4)))
#' 
#' @importFrom stats model.matrix
ARIMA <- function(formula, stepwise = TRUE, greedy = TRUE, 
                  approximation = FALSE, ...){
  arima_model$new(!!enquo(formula), stepwise = stepwise, greedy = greedy, 
                  approximation = approximation, ...)
}

#' @export
fitted.ARIMA <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.ARIMA <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
augment.ARIMA <- function(x, ...){
  x$est
}

#' @export
glance.ARIMA <- function(x, ...){
  x$fit
}

#' @export
tidy.ARIMA <- function(x, ...){
  x$par
}

#' @importFrom stats formula residuals
#' @export
forecast.ARIMA <- function(object, new_data = NULL, specials = NULL, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced")
  }
  
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
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
  
  # Output forecasts
  construct_fc(fc$pred, fc$se, dist_normal(fc$pred, fc$se))
}

#' @export
model_sum.ARIMA <- function(x){
  model_sum(x$model)
}

#' @export
model_sum.Arima <- function(x){
  order <- x$arma[c(1, 6, 2, 3, 7, 4, 5)]
  m <- order[7]
  result <- paste("ARIMA(", order[1], ",", order[2], ",", order[3], ")", sep = "")
  if (m > 1 && sum(order[4:6]) > 0) {
    result <- paste(result, "(", order[4], ",", order[5], ",", order[6], ")[", m, "]", sep = "")
  }
  if(length(x$coef) > sum(order[c(1,3,4,6)])){
    result <- paste("LM w/", result, "errors")
  }
  result
}

# Adjusted from robjhyndman/tsfeatures
#' @importFrom stats stl var
seas_heuristic <- function(x, period){
  stlfit <- stl(x, s.window = 13)
  remainder <- stlfit$time.series[,"remainder"]
  seasonal <- stlfit$time.series[, "seasonal"]
  vare <- var(remainder, na.rm = TRUE)
  max(0, min(1, 1 - vare/var(remainder + seasonal, na.rm = TRUE)))
}

# Adjusted from urca
#' @importFrom stats lm na.omit
kpss_test <- function(y, type = c("mu", "tau"), 
                       lags = c("short", "long", "nil"), use.lag = NULL){
  y <- na.omit(as.vector(y))
  n <- length(y)
  type <- match.arg(type)
  lags <- match.arg(lags)
  if (!(is.null(use.lag))) {
    lmax <- as.integer(use.lag)
    if (lmax < 0) {
      warning("\nuse.lag has to be positive and integer; lags='short' used.")
      lmax <- trunc(4 * (n/100)^0.25)
    }
  }
  else if (lags == "short") {
    lmax <- trunc(4 * (n/100)^0.25)
  }
  else if (lags == "long") {
    lmax <- trunc(12 * (n/100)^0.25)
  }
  else if (lags == "nil") {
    lmax <- 0
  }
  if (type == "mu") {
    cval <- as.matrix(t(c(0.347, 0.463, 0.574, 0.739)))
    colnames(cval) <- c("10pct", "5pct", "2.5pct", "1pct")
    rownames(cval) <- "critical values"
    res <- y - mean(y)
  }
  else if (type == "tau") {
    cval <- as.matrix(t(c(0.119, 0.146, 0.176, 0.216)))
    colnames(cval) <- c("10pct", "5pct", "2.5pct", "1pct")
    rownames(cval) <- "critical values"
    trend <- 1:n
    res <- residuals(lm(y ~ trend))
  }
  S <- cumsum(res)
  nominator <- sum(S^2)/n^2
  s2 <- sum(res^2)/n
  if (lmax == 0) {
    denominator <- s2
  }
  else {
    index <- 1:lmax
    x.cov <- sapply(index, function(x) t(res[-c(1:x)]) %*% 
                      res[-c((n - x + 1):n)])
    bartlett <- 1 - index/(lmax + 1)
    denominator <- s2 + 2/n * t(bartlett) %*% x.cov
  }
  teststat <- nominator/denominator
  list(y = y, type = type, lag = as.integer(lmax), 
      teststat = as.numeric(teststat), cval = cval)
}