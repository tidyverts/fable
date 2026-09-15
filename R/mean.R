#' @importFrom stats sd
train_mean <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    cli::cli_abort("Only univariate responses are supported by MEAN.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    cli::cli_abort("All observations are missing, a model cannot be estimated without data.")
  }

  n <- length(y)
  window_size <- specials$window[[1]]
  if (is.null(window_size)) {
    y_mean <- mean(y, na.rm = TRUE)
    fits <- rep(y_mean, n)
  }
  else {
    fits <- slide_dbl(y, mean,
      na.rm = TRUE,
      .size = window_size, .partial = TRUE
    )
    y_mean <- fits[length(fits)]
    fits <- dplyr::lag(fits)
  }
  res <- y - fits
  sigma <- sd(res, na.rm = TRUE)

  structure(
    list(
      fitted = fits,
      resid = res,
      mean = y_mean,
      sigma = sigma,
      nobs = sum(!is.na(y)),
      window = window_size %||% NA,
      time = list(start = unclass(.data)[[index_var(.data)]][[1]], interval = interval(.data))
    ),
    class = "model_mean"
  )
      
      # est = tibble(.fitted = fits, .resid = res),
      # fit = tibble(sigma2 = sigma^2),
      # spec = tibble(window_size = window_size %||% NA)
}

specials_mean <- new_specials(
  window = function(size = NULL) {
    size
  },
  .required_specials = "window"
)

#' Mean models
#'
#' \code{MEAN()} returns an iid model applied to the formula's response variable.
#'
#' @aliases report.model_mean
#'
#' @param formula Model specification.
#' @param ... Not used.
#'
#' @section Specials:
#'
#' \subsection{window}{
#' The `window` special is used to specify a rolling window for the mean.
#' \preformatted{
#' window(size = NULL)
#' }
#'
#' \tabular{ll}{
#'   `size`     \tab The size (number of observations) for the rolling window. If NULL (default), a rolling window will not be used.
#' }
#' }
#'
#' @return A model specification.
#'
#' @seealso
#' [Forecasting: Principles and Practices, Some simple forecasting methods (section 3.2)](https://otexts.com/fpp3/simple-methods.html)
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand))
#' @export
MEAN <- function(formula, ...) {
  mean_model <- new_model_class("mean",
    train = train_mean,
    specials = specials_mean
  )
  new_model_definition(mean_model, !!enquo(formula), ...)
}

#' @importFrom fabletools forecast
#' @importFrom stats qnorm time
#' @importFrom utils tail
#'
#' @inherit forecast.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   forecast()
#' @export
forecast.model_mean <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 5000, ...) {
  h <- NROW(new_data)

  y_mean <- object$mean
  n <- length(object$resid)
  sigma <- object$sigma

  # Produce forecasts
  if (bootstrap) { # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x) {
      generate(object, new_data, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose() %>%
      map(as.numeric)
    distributional::dist_sample(sim)
  } else {
    fc <- rep(y_mean, h)
    se <- sigma * sqrt(1 + 1 / n)
    distributional::dist_normal(fc, se)
  }
}

#' @inherit generate.ETS
#' @importFrom stats na.omit
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   generate()
#' @export
generate.model_mean <- function(x, new_data, bootstrap = FALSE, ...) {
  f <- x$mean

  if (!(".innov" %in% names(new_data))) {
    if (bootstrap) {
      res <- residuals(x)
      new_data$.innov <- sample(na.omit(res) - mean(res, na.rm = TRUE),
        NROW(new_data),
        replace = TRUE
      )
    }
    else {
      new_data$.innov <- stats::rnorm(NROW(new_data), sd = x$sigma)
    }
  }

  transmute(group_by_key(new_data), ".sim" := f + !!sym(".innov"))
}

#' @inherit interpolate.ARIMA
#'
#' @examples
#' library(tsibbledata)
#'
#' olympic_running %>%
#'   model(mean = MEAN(Time)) %>%
#'   interpolate(olympic_running)
#' @export
interpolate.model_mean <- function(object, new_data, specials, ...) {
  # Get inputs
  y <- new_data[[measured_vars(new_data)]]
  window_size <- object$window
  miss_val <- is.na(y)

  if (!is.na(window_size)) {
    fits <- dplyr::lag(
      slide_dbl(y, mean, na.rm = TRUE, .size = window_size, .partial = TRUE)
    )[miss_val]
  }
  else {
    fits <- object$mean
  }

  new_data[[measured_vars(new_data)]][miss_val] <- fits
  new_data
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   fitted()
#' @export
fitted.model_mean <- function(object, ...) {
  object$fitted
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   residuals()
#' @export
residuals.model_mean <- function(object, ...) {
  object$resid
}

#' Glance a average method model
#'
#' Construct a single row summary of the average method model.
#'
#' Contains the variance of residuals (`sigma2`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   glance()
#' @export
glance.model_mean <- function(x, ...) {
  tibble(sigma2 = x$sigma^2)
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' library(tsibbledata)
#' vic_elec %>%
#'   model(avg = MEAN(Demand)) %>%
#'   tidy()
#' @export
tidy.model_mean <- function(x, ...) {
  mu <- x$mean
  se <- x$sigma / sqrt(x$nobs)
  stat <- mu/se
  tibble(term = "mean", estimate = mu, std.error = se,
         statistic = stat,
         p.value = 2 * stats::pt(abs(stat), x$nobs - 1, lower.tail = FALSE))
}

#' @export
report.model_mean <- function(object, ...) {
  cat("\n")
  cat(paste("Mean:", round(object$mean, 4), "\n"))
  cat(paste("sigma^2:", round(object$sigma^2, 4), "\n"))
}

#' @export
model_sum.model_mean <- function(x) {
  paste0("MEAN") # , ", intToUtf8(0x3BC), "=", format(x$par$estimate))
}

#' Refit a MEAN model
#'
#' Applies a fitted average method model to a new dataset.
#'
#' @inheritParams refit.ARIMA
#' @param reestimate If `TRUE`, the mean for the fitted model will be re-estimated 
#' to suit the new data. 
#' 
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(MEAN(value))
#'
#' report(fit)
#'
#' fit %>%
#'   refit(lung_deaths_female) %>%
#'   report()
#' @export
refit.model_mean <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  # Update data for re-evaluation
  # update specials
  specials$window <- if(is.na(object$window)) NULL else object$window 

  if (reestimate) {
    return(train_mean(new_data, specials, ...))
  }
  
  y <- unclass(new_data)[[measured_vars(new_data)]]
  
  if (all(is.na(y))) {
    cli::cli_abort("All new observations are missing, model cannot be applied.")
  }

  if (!is_null(specials$window)) cli::cli_warn("A rolling mean model cannot be refitted, the most recent mean from the fitted model will be used as a fixed estimate of the mean.")
  
  n <- length(y)

  fits <- rep(object$mean, n)
  res <- y - fits
  sigma <- sd(res, na.rm = TRUE)
  
  object$fitted <- fits
  object$resid <- res
  object$sigma <- sigma
  object$nobs <- sum(!is.na(y))
  object$time <- list(start = unclass(new_data)[[index_var(new_data)]][[1]], interval = interval(new_data))
  object
}

#' Extend a fitted mean model with new data
#'
#' Applies a fitted average method model to a new (immediately subsequent)
#' portion of data, updating the model's fitted values, residuals and `sigma2`
#' without re-estimating the mean. Unlike [`refit.model_mean()`], `stream()`
#' does not need to reprocess the entire history, only the newly provided
#' observations, making it well suited to incrementally updating a model as
#' new observations arrive.
#'
#' @inheritParams refit.model_mean
#'
#' @details
#' For a fixed mean model, the `.fitted`, `.resid` and `sigma2` values
#' returned by `stream()` are identical to those obtained by refitting the
#' model (with the same fixed mean) to the complete series.
#'
#' For a rolling window mean model (specified with the `window()` special),
#' the rolling mean is continued across the new observations. The fitted values
#' and residuals are identical to those obtained by estimating the model on the
#' complete series, and the most recent rolling mean is retained for
#' forecasting.
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   filter(index < yearmonth("1979 Jan")) %>%
#'   model(MEAN(value))
#'
#' fit %>%
#'   stream(lung_deaths_male %>% filter(index >= yearmonth("1979 Jan"))) %>%
#'   report()
#' @export
stream.model_mean <- function(object, new_data, specials = NULL, ...) {
  n <- length(object$fitted)
  stream_start <- object$time$start + n * default_time_units(object$time$interval)
  if (unclass(new_data)[[index_var(new_data)]][1] != stream_start) {
    cli::cli_abort("Streaming to a mean model must start one step beyond the end of the trained data.")
  }

  y <- unclass(new_data)[[measured_vars(new_data)]]
  h <- length(y)
  window_size <- object$window

  if (is.na(window_size)) {
    fits <- rep(object$mean, h)
  } else {
    # Recover the last window of training observations (first fitted value is NA)
    y_past <- object$fitted + object$resid
    if (n > 1) y_past[1] <- object$fitted[2]
    y_past <- utils::tail(y_past, window_size)
    k <- length(y_past)

    # Continue the lagged rolling mean across the boundary
    roll <- slide_dbl(c(y_past, y), mean,
      na.rm = TRUE,
      .size = window_size, .partial = TRUE
    )
    fits <- roll[k + seq_len(h) - 1]
    object$mean <- roll[k + h]
  }
  res <- y - fits

  object$fitted <- c(object$fitted, fits)
  object$resid <- c(object$resid, res)
  object$sigma <- sd(object$resid, na.rm = TRUE)
  object$nobs <- object$nobs + sum(!is.na(y))

  object
}
