#' Random walk models
#' 
#' \code{RW()} returns a random walk model, which is equivalent to an ARIMA(0,1,0)
#' model with an optional drift coefficient included using \code{drift()}. \code{naive()} is simply a wrapper
#' to \code{rwf()} for simplicity. \code{snaive()} returns forecasts and
#' prediction intervals from an ARIMA(0,0,0)(0,1,0)m model where m is the
#' seasonal period.
#'
#' The random walk with drift model is \deqn{Y_t=c + Y_{t-1} + Z_t}{Y[t]=c +
#' Y[t-1] + Z[t]} where \eqn{Z_t}{Z[t]} is a normal iid error. Forecasts are
#' given by \deqn{Y_n(h)=ch+Y_n}{Y[n+h]=ch+Y[n]}. If there is no drift (as in
#' \code{naive}), the drift parameter c=0. Forecast standard errors allow for
#' uncertainty in estimating the drift parameter (unlike the corresponding
#' forecasts obtained by fitting an ARIMA model directly).
#'
#' The seasonal naive model is \deqn{Y_t= Y_{t-m} + Z_t}{Y[t]=Y[t-m] + Z[t]}
#' where \eqn{Z_t}{Z[t]} is a normal iid error.
#' 
#' @param data A data frame
#' @param formula Model specification.
#' 
#' @examples 
#' library(tsibbledata)
#' elecdemand %>% RW(Demand ~ drift())
#' 
#' @export
RW <- function(data, formula = ~ lag(1)){
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
  origin <- min(data[[expr_text(index(data))]])
  specials <- new_specials_env(
    lag = function(lag = 1){
      get_frequencies(lag, .data)
    },
    drift = function(drift = TRUE){
      drift
    },
    xreg = no_xreg,
    .env = caller_env(),
    .required_specials = c("lag"),
    .vals = list(.data = data, origin = origin)
  )
  
  estimate_RW(data = data, formula = formula, specials = specials, cl = cl)
}

#' @rdname RW
#'
#' @examples
#' 
#' Nile %>% as_tsibble %>% NAIVE
#'
#' @export
NAIVE <- RW

#' @rdname RW
#'
#' @examples
#' library(tsibbledata)
#' elecdemand %>% SNAIVE(Temperature ~ lag("day"))
#'
#' @export
SNAIVE <- function(data, formula = ~ lag("smallest")){
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
  origin <- min(data[[expr_text(index(data))]])
  specials <- new_specials_env(
    lag = function(lag = "smallest"){
      lag <- get_frequencies(lag, .data)
      if(lag == 1){
        abort("Non-seasonal model specification provided, use RW() or provide a different lag specification.")
      }
      lag
    },
    drift = function(drift = TRUE){
      drift
    },
    xreg = no_xreg,
    
    .env = caller_env(),
    .required_specials = c("lag"),
    .vals = list(.data = data, origin = origin)
  )
  
  estimate_RW(data=data, formula = formula, specials = specials, cl = cl)
}

#' @importFrom stats arima
estimate_RW <- function(data, formula, specials, cl){
  # Parse model
  model_inputs <- parse_model(data, formula, specials = specials)
  
  y <- eval_tidy(model_lhs(model_inputs$model), data = data)
  
  drift <- model_inputs$specials$drift[[1]] %||% FALSE
  lag <- model_inputs$specials$lag[[1]]
  
  y_na <- which(is.na(y))
  y_na <- y_na[y_na>lag]
  fits <- stats::lag(y, -lag)
  for(i in y_na){
    if(is.na(fits)[i]){
      fits[i] <- fits[i-lag]
    }
  }
  
  if(lag == 1 & !drift){
    method <- "NAIVE"
  }
  else if(lag != 1){
    method <- "SNAIVE"
  }
  else{
    method <- "RW"
  }
  
  fitted <- c(rep(NA, lag), utils::head(fits, -lag))
  if(drift){
    fit <- summary(stats::lm(y-fitted ~ 1, na.action=stats::na.exclude))
    b <- fit$coefficients[1,1]
    b.se <- fit$coefficients[1,2]
    sigma <- fit$sigma
    fitted <- fitted + b
    method <- paste(method, "w/ drift")
  }
  else{
    b <- b.se <- 0
    sigma <- stats::sd(y-fitted, na.rm=TRUE)
  }
  res <- y - fitted
  
  fit <- structure(
    list(
      par = tibble(term = "b", estimate = b, std.error = b.se),
      est = data %>% 
        transmute(
          !!model_lhs(model_inputs$model),
          .fitted = fitted,
          .resid = res
        ),
      fit = tibble(method = method,
                   lag = lag,
                   drift = drift,
                   sigma = sigma),
      future = mutate(new_data(data, lag), 
                      !!expr_text(model_lhs(model_inputs$model)) := utils::tail(fits, lag))
    ),
    class = "RW"
  )
  
  # Output model
  mable(
    data,
    model = fit,
    model_inputs
  )
}


#' @importFrom fablelite forecast
#' @importFrom stats qnorm time
#' @importFrom utils tail
#' @export
forecast.RW <- function(object, new_data = NULL, bootstrap = FALSE, times = 5000, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced")
  }
  
  h <- NROW(new_data)
  lag <- object$fit$lag
  fullperiods <- (h-1)/lag+1
  steps <- rep(1:fullperiods, rep(lag,fullperiods))[1:h]
  
  # Point forecasts
  fc <- rep(object$future[[measured_vars(object$future)[1]]], fullperiods)[1:h] +
    steps*object$par$estimate[1]
  
  # Intervals
  if (bootstrap){ # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x){
      simulate(object, new_data, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    dist <- sample_quantile(sim)
  }  else {
    mse <- mean(object$est$.resid^2, na.rm=TRUE)
    se  <- sqrt(mse*steps + (steps*object$par$std.error[1])^2)
    # Adjust prediction intervals to allow for drift coefficient standard error
    if (object$fit$drift) {
      se <- sqrt(se^2 + (seq(h) * object$par$std.error[1])^2)
    }
    dist <- new_fcdist(qnorm, fc, sd = se, abbr = "N")
  }
  
  construct_fc(new_data, fc, se, dist,
               expr_text(response(object)))
}


#' @export
simulate.RW <- function(object, new_data, bootstrap = FALSE, ...){
  if(!is_regular(new_data)){
    abort("Simulation new_data must be regularly spaced")
  }
  
  lag <- object$fit$lag
  fits <- select(rbind(object$est, object$future), !!index(object$est), !!response(object))
  start_idx <- min(new_data[[expr_text(index(new_data))]])
  start_pos <- match(start_idx, fits[[index(object$est)]])
  
  future <- fits[[expr_text(response(object))]][start_pos + seq_len(lag) - 1]
  
  if(any(is.na(future))){
    abort("The first lag window for simulation must be within the model's training set.")
  }
  
  if(is.null(new_data[[".innov"]])){
    if(bootstrap){
      new_data[[".innov"]] <- sample(stats::na.omit(object$est$.resid - mean(object$est$.resid, na.rm = TRUE)),
                                     NROW(new_data), replace = TRUE)
    }
    else{
      new_data[[".innov"]] <- stats::rnorm(NROW(new_data), sd = object$fit$sigma)
    }
  }

  sim_rw <- function(e){
    # Cumulate errors
    lag_grp <- rep_len(seq_len(lag), length(e))
    e <- split(e, lag_grp)
    cumulative_e <- unsplit(lapply(e, cumsum), lag_grp)
    future + cumulative_e 
  }
  
  new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := sim_rw(!!sym(".innov")))
}

#' @export
fitted.RW <- function(object, ...){
  select(object$est, ".fitted")
}

#' @export
residuals.RW <- function(object, ...){
  select(object$est, ".resid")
}

#' @export
augment.RW <- function(x, ...){
  x$est
}

#' @export
glance.RW <- function(x, ...){
  x$fit
}

#' @export
tidy.RW <- function(x, ...){
  x$par
}

#' @importFrom stats coef
#' @export
model_sum.RW <- function(x){
  x$fit$method
}