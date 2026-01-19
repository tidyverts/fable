#' Estimate an ARFIMA model
#'
#' Searches through the model space specified in the specials to identify a
#' suitable ARFIMA model. ARFIMA (AutoRegressive Fractionally Integrated Moving 
#' Average) models extend ARIMA models by allowing fractional differencing, 
#' which is useful for modeling long memory processes. The model is implemented
#' using [`fracdiff::fracdiff()`] and allows ARFIMA models to be used in the 
#' fable framework.
#'
#' @aliases report.fbl_ARFIMA
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param selection_metric A function used to compute a metric from the fitted
#' object which is minimised to select the best model.
#' @param stepwise,greedy,order_constraint,trace Arguments kept for API 
#' compatibility with `ARIMA()`. Currently not fully implemented for ARFIMA.
#' @param ... Further arguments passed to [`fracdiff::fracdiff()`].
#'
#' @section Parameterisation:
#'
#' An ARFIMA(p,d,q) model is defined as:
#'
#' \deqn{(1-\phi_1B - \cdots - \phi_p B^p)(1-B)^d (y_t - \mu) = (1 + \theta_1 B + \cdots + \theta_q B^q)\varepsilon_t}
#'
#' where \eqn{\mu} is the mean of the series, and \eqn{d} can take fractional 
#' values (typically between -0.5 and 0.5), allowing the model to capture long 
#' memory behavior. When \eqn{d} is an integer, the model reduces to a standard 
#' ARIMA model.
#' 
#' **Note:** This uses a mean form parameterisation where the data is de-meaned 
#' before fitting. This differs from [`ARIMA()`] which uses a constant form 
#' parameterisation.
#'
#' The fractional differencing operator \eqn{(1-B)^d} is computed using the 
#' fast algorithm of Jensen and Nielsen (2014), which is implemented in the 
#' fracdiff package.
#'
#' @section Specials:
#' 
#' The _specials_ define the space over which `ARFIMA` will search for the 
#' model that best fits the data. If the RHS of `formula` is left blank, the 
#' default search space is given by `pdq()`: a model with candidate 
#' non-seasonal terms and fractional differencing, but no exogenous regressors.
#'
#' Note that ARFIMA does not support seasonal differencing (PDQ terms). For 
#' seasonal data, consider using [`ARIMA()`] instead, or pre-process your data 
#' to remove seasonality.
#'
#' \subsection{pdq}{
#' The `pdq` special is used to specify the components of the ARFIMA model.
#' \preformatted{
#' pdq(p = 0:5, d = NULL, q = 0:5, 
#'     d_range = c(0, 0.5),
#'     p_init = 2, q_init = 2, fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `p`       \tab The order of the auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `d`       \tab The fractional differencing parameter. If `NULL` (default), it will be estimated. If a single numeric value is provided, it will be held fixed at that value. Unlike ARIMA, only a single value or `NULL` is allowed. \cr
#'   `q`       \tab The order of the moving average (MA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `d_range` \tab A numeric vector of length 2 specifying the range for estimating `d`. Only used when `d = NULL`. Typical values are between -0.5 and 0.5. \cr
#'   `p_init`  \tab If `stepwise = TRUE`, `p_init` provides the initial value for `p` for the stepwise search procedure. \cr
#'   `q_init`  \tab If `stepwise = TRUE`, `q_init` provides the initial value for `q` for the stepwise search procedure. \cr
#'   `fixed`   \tab A named list of fixed parameters for coefficients. The names identify the coefficient, beginning with either `ar` or `ma`, followed by the lag order. For example, `fixed = list(ar1 = 0.3, ma2 = 0)`.
#' }
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an ARFIMA model without explicitly 
#' using the `xreg()` special. Common exogenous regressor specials as specified 
#' in [`common_xregs`] can also be used. These regressors are handled using 
#' [stats::model.frame()], and so interactions and other functionality behaves 
#' similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows similar rules to 
#' [`stats::lm()`], where including `1` will add a constant and `0` or `-1` 
#' will remove the constant. If left out, the inclusion of a constant will be 
#' determined by minimising `ic`.
#'
#' \preformatted{
#' xreg(..., fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `...`   \tab Bare expressions for the exogenous regressors (such as `log(x)`) \cr
#'   `fixed` \tab A named list of fixed parameters for coefficients. The names identify the coefficient, and should match the name of the regressor. For example, `fixed = list(constant = 20)`.
#' }
#' }
#'
#' @seealso
#' [`ARIMA()`] for standard ARIMA models with integer differencing.
#' 
#' [Forecasting: Principles and Practices, ARIMA models (chapter 9)](https://otexts.com/fpp3/arima.html)
#' 
#' [`fracdiff::fracdiff()`] for the underlying fitting function.
#'
#' @references
#' Jensen, A. N. and Nielsen, M. Ø. (2014) A Fast Fractional Difference 
#' Algorithm. Journal of Time Series Analysis 35(5), 428–436. 
#' \doi{10.1111/jtsa.12074}
#'
#' @return A model specification.
#'
#' @examplesIf requireNamespace("fracdiff", quietly = TRUE)
#' library(tsibble)
#' library(dplyr)
#'
#' # Automatic ARFIMA specification
#' as_tsibble(sunspot.year) %>%
#'  model(arfima = ARFIMA(value)) %>%
#'  report()
#'  
#' # Manual ARFIMA specification with fixed d
#' as_tsibble(sunspot.year) %>%
#'   model(arfima = ARFIMA(value ~ pdq(p = 1, d = 0.3, q = 1))) %>%
#'   report()
#' 
#' @export
ARFIMA <- function(
  formula,
  ic = c("aicc", "aic", "bic"),
  selection_metric = function(x) x[[ic]],
  stepwise = TRUE,
  greedy = TRUE,
  order_constraint = p + q <= 6,
  trace = FALSE,
  ...
) {
  ic <- match.arg(ic)
  stopifnot(is.function(selection_metric))

  arfima_model <- new_model_class(
    "ARFIMA",
    train = train_arfima,
    specials = specials_arfima,
    origin = NULL,
    check = all_tsbl_checks
  )

  new_model_definition(
    arfima_model,
    !!enquo(formula),
    ic = ic,
    selection_metric = selection_metric,
    stepwise = stepwise,
    greedy = greedy,
    order_constraint = enexpr(order_constraint),
    trace = trace,
    ...
  )
}

specials_arfima <- new_specials(
  pdq = function(
    p = 0:5,
    d = NULL,
    q = 0:5,
    d_range = c(0, 0.5),
    p_init = 2,
    q_init = 2,
    fixed = list()
  ) {
    if (self$stage %in% c("estimate", "refit")) {
      p <- p[p <= floor(NROW(self$data) / 3)]
      q <- q[q <= floor(NROW(self$data) / 3)]
    }

    p_init <- p[which.min(abs(p - p_init))]
    q_init <- q[which.min(abs(q - q_init))]

    if (!is.null(d) && length(d) != 1 && !is.numeric(d)) {
      abort("For ARFIMA, `d` (if supplied) must be a single numeric value.")
    }

    # Limit drange if d is specified
    if (!is.null(d)) {
      d_range <- c(d, d)
    }

    if (!all(grepl("^(ma|ar)\\d+", names(fixed)))) {
      abort(
        "The 'fixed' coefficients for pdq() must begin with ar or ma, followed by a lag number."
      )
    }

    as.list(environment())
  },

  # ARFIMA has no seasonal PDQ special.
  PDQ = function(...) {
    abort("Seasonal PDQ() is not currently supported for ARFIMA models.")
  },

  # xreg handling can be the same as for ARIMA
  common_xregs,

  xreg = function(..., fixed = list()) {
    dots <- enexprs(...)
    env <- map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
    env <- if (!is_empty(env)) get_env(env[[1]]) else base_env()
    constants <- map_lgl(dots, inherits, "numeric")
    constant_forced <- any(map_lgl(dots[constants], `%in%`, 1))
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    env <- env_bury(env, lag = lag)

    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    tm <- terms(xreg)
    constant <- as.logical(tm %@% "intercept")
    xreg <- model.matrix(tm, xreg)
    if (constant) {
      xreg <- xreg[, -1, drop = FALSE]
    }

    list(
      constant = if (!constant || constant_forced) constant else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg,
      fixed = fixed
    )
  },

  .required_specials = c("pdq"),
  .xreg_specials = names(common_xregs)
)

# ' Fractional differences
# ' 
# ' Computes the fractional difference of a numeric vector using the fast algorithm by Jensen and Nielsen (2014).
# ' 
# ' @param x Numeric vector to difference.
# ' @param d Fractional differencing parameter.
# ' 
# ' @return A numeric vector representing the fractional difference of `x`.
# ' 
# ' @references Jensen, Andreas Noack and Nielsen, Morten Ørregaard (2014) A Fast Fractional Difference Algorithm. Journal of Time Series Analysis 35(5), 428–436; [doi:10.1111/jtsa.12074](https://doi.org/10.1111/jtsa.12074).
# ' 
# ' @seealso [fracdiff::diffseries()] for the original function which additionally differences the data.
# ' 
#' @importFrom stats fft nextn
fracdiff <- function(x, d) {
  stopifnot((iT <- length(x)) >= 2)
  np2 <- nextn(iT + iT - 1L)
  # Pad with zero
  pad <- rep.int(0, np2 - iT)
  k <- seq_len(iT - 1L)
  b <- c(1, cumprod((k - (d+1))/ k), pad)
  ## ~= convolve(x, b, type = "filter") :
  dx <- fft(fft(b) * fft(c(x, pad)), inverse =TRUE)[seq_len(iT)] / np2
  Re(dx)
}

# ' Fractional Integration: Inverse of Fractional Differencing
# ' 
# ' Computes the inverse function of the fractional differencing function [fracdiff()].
# '
# ' @param x Numeric vector to integrate.
# ' @param d Fractional integration parameter.
# ' @param xi A numeric vector containing the initial values for the fractional integral.
# ' 
# ' @return A numeric vector representing the fractional integral of `x`.
# ' 
# ' @examples
# ' 
# ' s <- 1:10
# ' d <- fracdiff(s, 0.5)
# ' fracdiffinv(d, 0.5)
fracdiffinv <- function (x, d, xi = NULL) {
  # Trim xi to relevant history based on FFT padding used in fracdiff
  # The FFT uses np2 = nextn(2*n - 1), so we only need approximately n values
  if (!is.null(xi) && length(xi) > length(x)) {
    xi <- tail(xi, length(x))
  }
  
  fracdiff(c(xi, x), -d)[seq(length(xi) + 1L, length.out = length(x))]
}

#' @importFrom stats ts
train_arfima <- function(
  .data,
  specials,
  ... # Passed through to `train_arima()`
) {
  # Requires fracdiff for ARFIMA fitting
  require_package("fracdiff")

  # Extract response
  y <- unclass(.data)[[measured_vars(.data)]]
  idx <- unclass(.data)[[index_var(.data)]]
  
  # Demean data
  y_mu <- mean(y, na.rm = TRUE)
  y <- y - y_mu
  
  # Initialise fractional differencing `d` with `ARFIMA(2, d, 0)` fit
  if (is.null(specials$pdq[[1]]$d)) {
    d_init <- fracdiff::fracdiff(y, nar = 2, nma = 0)$d
  } else 
  {
    d_init <- specials$pdq[[1L]]$d
  }
  
  # Fractionally difference the data
  yd <- fracdiff(y, d_init)
  
  # Select ARFIMA hyperparameters p,q with auto.arima
  specials$pdq[[1L]]$d <- 0L
  specials$PDQ <- list(list(period = 1L, P = 0L, D = 0L, Q = 0L)) # No seasonal part for ARFIMA
  specials$xreg <- specials$xreg %||% list(list(constant = FALSE)) 
  .data[[measured_vars(.data)]] <- yd
  fit <- train_arima(.data, specials, ...)

  # Refit ARFIMA with selected p,q
  # TODO: handle xregs and constant before fracdiff::fracdiff?
  p <- specials$pdq[[1L]]$p <- fit$spec$p
  q <- specials$pdq[[1L]]$q <- fit$spec$q
  c <- specials$xreg[[1L]]$constant <- fit$spec$constant
  fit <- fracdiff::fracdiff(
    y, nar = p, nma = q,
    drange = specials$pdq[[1L]]$d_range
  )
  d <- fit$d

  # Refine ARIMA coefficients with MLE
  yd <- fracdiff(y, fit$d)
  .data[[measured_vars(.data)]] <- yd
  fit <- train_arima(.data, specials, ...)
  fit$spec$d <- d
  fit$spec$mu <- y_mu

  class(fit) <- c("fbl_ARFIMA", class(fit))
  fit
}

#' @export
model_sum.fbl_ARFIMA <- function(x) {
  sprintf(
    "ARFIMA(%i,%.2f,%i)%s",
    x$spec$p,
    x$spec$d,
    x$spec$q,
    if (isTRUE(x$spec$constant)) " w/ mean" else ""
  )
}

#' @export
fitted.fbl_ARFIMA <- function(object, ...) {

  # @RH: The {forecast} package does strange things with the fitted residuals
  #      (perhaps it doesn't fractionally integrate them?)
  fracdiffinv(
    NextMethod(),
    d = object$spec$d,
    # xi = ??? -- are there suitable non-zero initial values for this?
  ) + object$spec$mu
}

#' @importFrom stats frequency
#' @export
forecast.fbl_ARFIMA <- function(
  object,
  new_data = NULL,
  specials = NULL,
  bootstrap = FALSE,
  times = 5000,
  ...
) {
  # Restore training data to initialise fractional integration
  y_train <- object$est$.fitted + object$est$.resid
  n <- length(y_train)

  # Obtain ARMA forecasts
  fc <- NextMethod()
  h <- length(fc)
  
  # Extract ARMA coefficients
  p <- object$spec$p
  q <- object$spec$q
  phi <- theta <- numeric(h)
  phi[seq_len(p)] <- object$model$coef[grep("^ar", names(object$model$coef))]
  theta[seq_len(q)] <- object$model$coef[grep("^ma", names(object$model$coef))]
  
  # Binomial coefficient for expansion of d
  bin.c <- (-1)^(0:(n + h)) * choose(object$spec$d, (0:(n + h)))
  
  # Calculate psi weights
  new.phi <- psi <- numeric(h)
  psi[1] <- new.phi[1] <- 1
  if (h > 1) {
    new.phi[2:h] <- -bin.c[2:h]
    for (i in 2:h) {
      if (p > 0) {
        new.phi[i] <- sum(phi[1:(i - 1)] * bin.c[(i - 1):1]) - bin.c[i]
      }
      psi[i] <- sum(new.phi[2:i] * rev(psi[1:(i - 1)])) + theta[i - 1]
    }
  }

  distributional::dist_normal(
    # Fractionally integrated mean
    mu = fracdiffinv(mean(fc), d = object$spec$d, xi = y_train) + object$spec$mu,
    # Approximate standard error via psi-weights
    sigma = sqrt(cumsum(psi^2) * object$fit$sigma2)
  )
}

# A very simple generate method: simulate via forecast.fracdiff with simulate=TRUE or bootstrap.
#' @export
generate.fbl_ARFIMA <- function(
  x,
  new_data,
  specials,
  bootstrap = FALSE,
  times = 1,
  ...
) {
  # Simulate from ARMA model
  .sim <- NextMethod()

  # Initial values based on `new_data` time points
  y_init <- x$est$.fitted + x$est$.resid

  # ASSUME: simulation time points are contiguous from minimum time
  t_init <- min(new_data[[index_var(new_data)]])
  
  # Truncate initial values to those before t_init
  t_chronon <- sum(unlist(x$tsp$interval))
  y_init <- y_init[seq_len((t_init - x$tsp$range[1])/t_chronon)]

  # Fractionally integrate simulated data
  transmute(
    group_by_key(.sim), 
    ".sim" := fracdiffinv(x = .data$.sim, d = x$spec$d, xi = y_init) + x$spec$mu
  )
}

#' @export
refit.fbl_ARFIMA <- function(
  object,
  new_data,
  specials = NULL,
  reestimate = FALSE,
  ...
) {
  # Re-estimate coefficients using train_arfima()
  if (reestimate) {
    require_package("fracdiff")
    specials$pdq[[1]][c("p", "d", "q", "p_init", "q_init")] <- 
      as.list(object$spec[c("p", "d", "q", "p", "q")])
    
    # Re-estimate ARFIMA coefficients
    return(
      train_arfima(new_data, specials, ...)
    )
  }

  # Fractionally difference the new data
  y <- unclass(new_data)[[measured_vars(new_data)]]
  yd <- fracdiff(y - object$spec$mu, object$spec$d)
  new_data[[measured_vars(new_data)]] <- yd

  # Refit ARMA model
  fit <- NextMethod()
  
  # Class the results
  fit$spec$d <- object$spec$d
  fit$spec$mu <- object$spec$mu
  class(fit) <- c("fbl_ARFIMA", class(fit))
  fit
}
