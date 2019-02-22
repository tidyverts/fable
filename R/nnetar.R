#' @importFrom stats ar complete.cases
train_nnetar <- function(.data, formula, specials, size, repeats, scale_inputs, ...){
  require_package("nnet")
  
  y <- x <- .data[[measured_vars(.data)]]
  n <- length(x)
  
  if (n < 3) {
    stop("Not enough data to fit a model")
  }
  
  # Get args
  p <- P <- period <- NULL 
  assignSpecials(specials["AR"])
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  
  # Check for constant data in time series
  constant_data <- is.constant(x)
  if (constant_data){
    warn("Constant data, setting p=1, P=0, scale.inputs=FALSE")
    scale.inputs <- FALSE
    p <- 1
    P <- 0
  }
  ## Check for constant data in xreg
  if (!is.null(xreg)){
    xreg <- as.matrix(xreg)
    if (any(apply(xreg, 2, is.constant))){
      warning("Constant xreg column, setting scale_inputs=FALSE")
      scale_inputs <- FALSE
    }
  }
  
  # Scale inputs
  y_scale <- xreg_scale <- NULL
  if(scale_inputs){
    x <- scale(x, center = TRUE, scale = TRUE)
    # y_scale <- tibble(
    #   term = c("y_center", "y_scale"),
    #   estimate = c(attr(x, "scaled:center"), attr(x, "scaled:scale"))
    # )
    y_scale <- list(
      center = attr(x, "scaled:center"), 
      scale = attr(x, "scaled:scale")
    )
    x <- c(x)
    
    if(!is.null(xreg)){
      xreg <- scale(xreg, center = TRUE, scale = TRUE)
      # xreg_scale <- tibble(
      #   term = c("xreg_center", "xreg_scale"),
      #   estimate = c(attr(xreg, "scaled:center"), attr(xreg, "scaled:scale"))
      # )
      xreg_scale <- list(
        center = attr(xreg, "scaled:center"), 
        scale = attr(xreg, "scaled:scale"))
    }
  }
  
  # Construct lagged matrix
  if(is.null(p)){
    if(period > 1 && n >= 2 * period){
      y_sa <- stl(ts(x, frequency = period), s.window = 13)
      y_sa <- y_sa$time.series[,"trend"] + y_sa$time.series[,"remainder"]
    }
    else{
      y_sa <- x
    }
    p <- max(length(ar(y_sa)$ar), 1)
  }
  if (p >= n) {
    warn("Reducing number of lagged inputs due to short series")
    p <- n - 1
  }
  lags <- 1:p
  if (P > 0 && n >= period * P + 2) {
    lags <- sort(unique(c(lags, period * (1:P))))
  }
  
  maxlag <- max(lags)
  nlag <- length(lags)
  x_lags <- matrix(NA_real_, ncol = nlag, nrow = n - maxlag)
  for (i in 1:nlag)
    x_lags[, i] <- x[(maxlag - lags[i] + 1):(n - lags[i])]
  x <- x[-(1:maxlag)]
  # Add xreg into lagged matrix
  x_lags <- cbind(x_lags, xreg[-(1:maxlag), ])
  if (is.null(size)) {
    size <- round((NCOL(x_lags) + 1) / 2)
  }
  
  # Remove missing values if present
  j <- complete.cases(x_lags, x)
  x_lags <- x_lags[j, , drop=FALSE]
  x <- x[j]
  ## Stop if there's no data to fit
  if (NROW(x_lags) == 0) {
    abort("No data to fit (possibly due to missing values)")
  }
  
  # Fit the nnet
  nn_models <- map(seq_len(repeats),
                   function(.){wrap_nnet(x_lags, x, size = size, ...)})
  
  # Construct model output
  fits <- map(nn_models, predict) %>% 
    transpose %>% 
    map_dbl(function(x) mean(unlist(x)))
  fits <- c(rep(NA_real_, maxlag), fits)
  if(scale_inputs){
    fits <- fits * y_scale$scale + y_scale$center
  }
  res <- y - fits
  
  structure(
    list(
      model = nn_models,
      par = tibble(),
      est = .data %>% 
        mutate(
          .fitted = fits,
          .resid = res
        ),
      fit = tibble(sigma = sd(res, na.rm = TRUE)),
      spec = tibble(period = period, p = p, P = P, size = size, lags = list(lags)),
      scales = list(y = y_scale, xreg = xreg_scale),
      future = mutate(new_data(.data, maxlag), 
                      !!expr_text(model_lhs(self)) := utils::tail(x, maxlag))
    ),
    class = "NNETAR"
  )
}

# Wrap nnet to change the default for linout to be TRUE
wrap_nnet <- function(x, y, repeats, linout = TRUE, trace = FALSE, ...) {
  nnet::nnet(x, y, linout = linout, trace = trace, ...)
}

specials_nnetar <- new_specials(
  AR = function(p = NULL, P = 1, period = "smallest"){
    period <- get_frequencies(period, self$data)
    if (period == 1) {
      if(!missing(P) && P > 0){
        warn("Non-seasonal data, ignoring seasonal lags")
      }
      P <- 0
    }
    if (P > 0 && NROW(self$data) < period * P + 2) {
      warn("Series too short for seasonal lags")
      P <- 0
    }
    as.list(environment())
  },
  common_xregs,
  xreg = model_xreg,
  .required_specials = c("AR")
)

nnetar_model <- R6::R6Class(NULL,
                        inherit = fablelite::model_definition,
                        public = list(
                          model = "NNETAR",
                          train = train_nnetar,
                          specials = specials_nnetar,
                          origin = NULL
                        )
)

#' Neural Network Time Series Forecasts
#'
#' Feed-forward neural networks with a single hidden layer and lagged inputs
#' for forecasting univariate time series.
#'
#' A feed-forward neural network is fitted with lagged values of the response as
#' inputs and a single hidden layer with `size` nodes. The inputs are for
#' lags 1 to `p`, and lags `m` to `mP` where
#' `m` is the seasonal period specified.
#' 
#' If exogenous regressors are provided, its columns are also used as inputs.
#' Missing values are currently not supported by this model.
#' A total of `repeats` networks are
#' fitted, each with random starting weights. These are then averaged when
#' computing forecasts. The network is trained for one-step forecasting.
#' Multi-step forecasts are computed recursively.
#'
#' For non-seasonal data, the fitted model is denoted as an NNAR(p,k) model,
#' where k is the number of hidden nodes. This is analogous to an AR(p) model
#' but with non-linear functions. For seasonal data, the fitted model is called
#' an NNAR(p,P,k)\[m\] model, which is analogous to an ARIMA(p,0,0)(P,0,0)\[m\]
#' model but with non-linear functions.
#'
#' @inheritParams MEAN
#' @param size Number of nodes in the hidden layer. Default is half of the
#' number of input nodes (including external regressors, if given) plus 1.
#' @param repeats Number of networks to fit with different random starting
#' weights. These are then averaged when producing forecasts.
#' @param scale_inputs If TRUE, inputs are scaled by subtracting the column
#' means and dividing by their respective standard deviations. Scaling is
#' applied after transformations.
#' @param ... Other arguments passed to `\link[nnet]{nnet}`.
#'
#' @author Gabriel Caceres, Mitchell O'Hara-Wild & Rob J Hyndman
#' 
#' @export
NNETAR <- function(formula, size = NULL, repeats = 20, scale_inputs = TRUE, ...){
  nnetar_model$new(!!enquo(formula), size = size, repeats = repeats,
                   scale_inputs = scale_inputs, ...)
}

#' @export
forecast.NNETAR <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 1000, ...){
  require_package("nnet")
  
  # Prepare xreg
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  if(!is.null(xreg)){
    xreg <- as.matrix(xreg)
    if(!is.null(object$scales$xreg)){
      xreg <- scale(xreg, center = object$scales$xreg$center, scale = object$scales$xreg$scale)
    }
  }
  
  # Extract model attributes
  lags <- object$spec$lags[[1]]
  maxlag <- max(lags)
  future_lags <- rev(object$future[[measured_vars(object$future)]])
  
  # Compute 1-step forecasts iteratively
  h <- NROW(new_data)
  fc <- numeric(h)
  for (i in 1:h)
  {
    fcdata <- c(future_lags[lags], xreg[i, ])
    if (any(is.na(fcdata))) {
      abort("I can't use NNETAR to forecast with missing values near the end of the series.")
    }
    fc[i] <- mean(map_dbl(object$model, predict, newdata = fcdata))
    future_lags <- c(fc[i], future_lags[-maxlag])
  }
  # Re-scale point forecasts
  if (!is.null(object$scales$y)) {
    fc <- fc * object$scales$y$scale + object$scales$y$center
  }
  
  # Compute forecast intervals
  sim <- map(seq_len(times), function(x){
    imitate(object, new_data, specials = specials, bootstrap = bootstrap)[[".sim"]]
  }) %>%
    transpose %>%
    map(as.numeric)
  se <- map_dbl(sim, stats::sd)
  dist <- dist_sim(sim)
  
  construct_fc(fc, se, dist)
}

#' @export
imitate.NNETAR <- function(object, new_data, specials = NULL, bootstrap = FALSE, ...){
  # Prepare xreg
  xreg <- specials[c("xreg", names(common_xregs))] %>% 
    compact() %>% 
    map(function(.x){invoke("cbind", .x)}) %>% 
    invoke("cbind", .)
  if(!is.null(xreg)){
    xreg <- as.matrix(xreg)
    if(!is.null(object$scales$xreg)){
      xreg <- scale(xreg, center = object$scales$xreg$center, scale = object$scales$xreg$scale)
    }
  }
  
  if(is.null(new_data[[".innov"]])){
    if(bootstrap){
      res <- stats::na.omit(object$est[[".resid"]] - mean(object$est[[".resid"]], na.rm = TRUE))
      if (!is.null(object$scales$y)) {
        res <- res / object$scales$y$scale
      }
      new_data[[".innov"]] <- sample(res, NROW(new_data), replace = TRUE)
    }
    else{
      if (!is.null(object$scales$y)) {
        sigma <- sd(object$est[[".resid"]] / object$scales$y$scale, na.rm = TRUE)
      }
      else{
        sigma <- object$fit$sigma
      }
      new_data[[".innov"]] <- stats::rnorm(NROW(new_data), sd = sigma)
    }
  }
  else{
    if (!is.null(object$scales$y)) {
      new_data[[".innov"]] <- new_data[[".innov"]] / object$scales$y$scale
    }
  }
  
  # Extract model attributes
  lags <- object$spec$lags[[1]]
  maxlag <- max(lags)
  future_lags <- rev(object$future[[measured_vars(object$future)]])
  
  sim_nnetar <- function(e){
    path <- numeric(length(e))
    for (i in 1:length(e))
    {
      fcdata <- c(future_lags[lags], xreg[i, ])
      if (any(is.na(fcdata))) {
        abort("I can't use NNETAR to forecast with missing values near the end of the series.")
      }
      path[i] <- mean(map_dbl(object$model, predict, newdata = fcdata)) + e[i]
      future_lags <- c(path[i], future_lags[-maxlag])
    }
    
    # Re-scale simulated paths
    if (!is.null(object$scales$y)) {
      path <- path * object$scales$y$scale + object$scales$y$center
    }
    path
  }

  new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := sim_nnetar(!!sym(".innov")))
}

#' @export
fitted.NNETAR <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.NNETAR <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
augment.NNETAR <- function(x, ...){
  x$est
}

#' @export
glance.NNETAR <- function(x, ...){
  x$fit
}

#' @export
tidy.NNETAR <- function(x, ...){
  x$par
}

#' @export
report.NNETAR <- function(object, ...){
  cat(paste("\nAverage of", length(object$model), "networks, each of which is\n"))
  print(object$model[[1]])
  cat(
    "\nsigma^2 estimated as ", 
    format(mean(residuals(object) ^ 2, na.rm = TRUE), digits = 4),
    "\n", sep = ""
  )
  invisible(object)
}

#' @importFrom stats coef
#' @export
model_sum.NNETAR <- function(x){
  p <- x$spec$p
  P <- x$spec$P
  sprintf("NNAR(%s,%i)%s",
          ifelse(P > 0, paste0(p, ",", P), p),
          x$spec$size,
          ifelse(P > 0, paste0("[", x$spec$period, "]"), "")
  )
}