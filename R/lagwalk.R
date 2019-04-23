train_lagwalk <- function(.data, formula, specials, ...){
  if(length(measured_vars(.data)) > 1){
    abort("Only univariate responses are supported by lagwalks.")
  }
  
  y <- .data[[measured_vars(.data)]]
  n <- length(y)
  
  if(all(is.na(y))){
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  drift <- specials$drift[[1]] %||% FALSE
  lag <- specials$lag[[1]]
  
  y_na <- which(is.na(y))
  y_na <- y_na[y_na>lag]
  fits <- stats::lag(y, -lag)
  for(i in y_na){
    if(is.na(fits)[i]){
      fits[i] <- fits[i-lag]
    }
  }
  
  fitted <- c(rep(NA, min(lag, n)), utils::head(fits, -lag))
  if(drift){
    fit <- summary(stats::lm(y-fitted ~ 1, na.action=stats::na.exclude))
    b <- fit$coefficients[1,1]
    b.se <- fit$coefficients[1,2]
    sigma <- fit$sigma
    fitted <- fitted + b
  }
  else{
    b <- b.se <- NULL
    sigma <- stats::sd(y-fitted, na.rm=TRUE)
  }
  res <- y - fitted
  
  structure(
    list(
      par = tibble(term = if(drift) "b" else chr(),
                   estimate = b%||%dbl(), std.error = b.se%||%dbl()),
      est = tibble(
          .fitted = fitted,
          .resid = res
        ),
      fit = tibble(sigma = sigma),
      spec = tibble(lag = lag, drift = drift),
      future = mutate(new_data(.data, lag), 
                      !!expr_text(model_lhs(self)) := y[c(rep(NA, max(0, lag - n)), y[seq_len(min(n, lag)) + n - min(n, lag)])]
      )
    ),
    class = "RW"
  )
}

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
#' @param formula Model specification (see "Specials" section).
#' @param ... Not used.
#' 
#' 
#' @section Specials:
#' 
#' \subsection{lag}{
#' The `lag` special is used to specify the lag order for the random walk process. 
#' If left out, this special will automatically be included.
#' 
#' \preformatted{
#' lag(lag = NULL)
#' }
#'
#' \tabular{ll}{
#'   `lag`        \tab The lag order for the random walk process. If `lag = m`, forecasts will return the observation from `m` time periods ago. This can also be provided as text indicating the duration of the lag window (for example, annual seasonal lags would be "1 year").
#' }
#' }
#'  
#' \subsection{drift}{
#' The `drift` special can be used to include a drift/trend component into the model. By default, drift is not included unless `drift()` is included in the formula.
#' 
#' \preformatted{
#' drift(drift = TRUE)
#' }
#' 
#' \tabular{ll}{
#'   `drift`      \tab If `drift = TRUE`, a drift term will be included in the model.
#' }
#' }
#' 
#' @examples 
#' library(tsibbledata)
#' aus_production %>% 
#'   model(rw = RW(Beer ~ drift()))
#' 
#' @export
RW <- function(formula, ...){
  rw_model <- new_model_class("RW", train = train_lagwalk,
                              specials = new_specials(
                                lag = function(lag = NULL){
                                  if(is.null(lag)){
                                    lag <- 1
                                  }
                                  get_frequencies(lag, self$data, .auto = "smallest")
                                },
                                drift = function(drift = TRUE){
                                  drift
                                },
                                xreg = no_xreg,
                                .required_specials = c("lag")
                              ),
                              check = all_tsbl_checks)
  new_model_definition(rw_model, !!enquo(formula), ...)
}

#' @rdname RW
#'
#' @examples
#' 
#' Nile %>% as_tsibble %>% model(NAIVE())
#'
#' @export
NAIVE <- RW

#' @rdname RW
#'
#' @examples
#' library(tsibbledata)
#' aus_production %>% 
#'   model(snaive = SNAIVE(Beer ~ lag("year")))
#'
#' @export
SNAIVE <- function(formula, ...){
  snaive_model <- new_model_class("RW", train = train_lagwalk,
                              specials = new_specials(
                                lag = function(lag = NULL){
                                  lag <- get_frequencies(lag, self$data, .auto = "smallest")
                                  if(lag == 1){
                                    abort("Non-seasonal model specification provided, use RW() or provide a different lag specification.")
                                  }
                                  lag
                                },
                                drift = function(drift = TRUE){
                                  drift
                                },
                                xreg = no_xreg,
                                .required_specials = c("lag")
                              ),
                              check = all_tsbl_checks)
  new_model_definition(snaive_model, !!enquo(formula), ...)
}

#' @importFrom fablelite forecast
#' @importFrom stats qnorm time
#' @importFrom utils tail
#' @export
forecast.RW <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 5000, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced")
  }
  
  h <- NROW(new_data)
  lag <- object$spec$lag
  fullperiods <- (h-1)/lag+1
  steps <- rep(1:fullperiods, rep(lag,fullperiods))[1:h]
  
  b <- object$par$estimate
  b.se <- object$par$std.error
  if(!object$spec$drift){
    b <- b.se <- 0
  }
  
  # Point forecasts
  fc <- rep(object$future[[measured_vars(object$future)[1]]], fullperiods)[1:h] +
    steps*b
  
  # Intervals
  if (bootstrap){ # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x){
      imitate(object, new_data, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    dist <- dist_sim(sim)
  }  else {
    mse <- mean(object$est$.resid^2, na.rm=TRUE)
    if(is.nan(mse)) mse <- NA
    se  <- sqrt(mse*steps + (steps*b.se)^2)
    # Adjust prediction intervals to allow for drift coefficient standard error
    if (object$spec$drift) {
      se <- sqrt(se^2 + (seq(h) * b.se)^2)
    }
    dist <- dist_normal(fc, se)
  }
  
  construct_fc(fc, se, dist)
}


#' @export
imitate.RW <- function(object, new_data, bootstrap = FALSE, ...){
  if(!is_regular(new_data)){
    abort("Simulation new_data must be regularly spaced")
  }
  
  lag <- object$spec$lag
  if(object$spec$drift){
    b <- object$par$estimate
  } else {
    b <- 0
  }
  fits <- select(rbind(object$est, object$future), !!index(object$est), !!measured_vars(object$future))
  start_idx <- min(new_data[[expr_text(index(new_data))]])
  start_pos <- match(start_idx, fits[[index(object$est)]])
  
  future <- fits[[measured_vars(object$future)]][start_pos + seq_len(lag) - 1]
  
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
    dx <- e + b
    lag_grp <- rep_len(seq_len(lag), length(dx))
    dx <- split(dx, lag_grp)
    cumulative_e <- unsplit(lapply(dx, cumsum), lag_grp)
    rep_len(future, length(dx)) + cumulative_e
  }
  
  new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := sim_rw(!!sym(".innov")))
}

#' @export
fitted.RW <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.RW <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
glance.RW <- function(x, ...){
  x$fit
}

#' @export
tidy.RW <- function(x, ...){
  x$par
}

#' @export
report.RW <- function(object, ...){
  cat("\n")
  if (object$spec$drift) {
    row <- object$par$term == "b"
    cat(paste("Drift: ", round(object$par$estimate[row], 4), 
              " (se: ", round(object$par$std.error[row], 4), ")\n", sep = ""))
  }
  cat(paste("Residual sd:", round(object$fit$sigma, 4), "\n"))
}

#' @importFrom stats coef
#' @export
model_sum.RW <- function(x){
  if(x$spec$lag == 1 & !x$spec$drift){
    method <- "NAIVE"
  }
  else if(x$spec$lag != 1){
    method <- "SNAIVE"
  }
  else{
    method <- "RW"
  }
  if(x$spec$drift){
    method <- paste(method, "w/ drift")
  }
  method
}