#' @importFrom stats sd
train_mean <- function(.data, specials, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by MEAN.")
  }

  y <- unclass(.data)[[measured_vars(.data)]]

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
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
      window = window_size %||% NA
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
#' The model does not use any specials, and so everything on the formula's
#' right-hand-side will be ignored.
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
    se <- map_dbl(sim, stats::sd)
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

  new_data %>%
    group_by_key() %>%
    transmute(".sim" := f + !!sym(".innov"))
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
