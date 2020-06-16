globalVariables(c("p", "P", "d", "D", "q", "Q", "constant"))

#' @importFrom stats approx lm ts
train_arima <- function(.data, specials, ic = "aicc",
                        stepwise = TRUE, greedy = TRUE, approximation = NULL,
                        order_constraint = p + q + P + Q <= 6 & (!constant | d + D < 2),
                        unitroot_spec = unitroot_options(),
                        fixed = NULL, ...) {
  if (length(measured_vars(.data)) > 1) {
    abort("Only univariate responses are supported by ARIMA.")
  }
  
  # Get args
  if(length(specials$pdq) > 1 || length(specials$PDQ) > 1){
    warn("Only one special for `pdq()` and `PDQ()` is allowed, defaulting to the first usage")
  }
  pdq <- specials$pdq[[1]]
  PDQ <- specials$PDQ[[1]]
  period <- PDQ$period
  special_fixed <- c(specials$pdq[[1]]$fixed, specials$PDQ[[1]]$fixed, specials$xreg[[1]]$fixed)

  # Get response
  y <- x <- ts(unclass(.data)[[measured_vars(.data)]], frequency = period)

  if (all(is.na(y))) {
    abort("All observations are missing, a model cannot be estimated without data.")
  }

  # Get xreg
  constant <- specials$xreg[[1]]$constant %||% c(TRUE, FALSE)
  xreg <- specials$xreg[[1]]$xreg

  # Check xreg
  if (!is_empty(xreg)) {
    # Now check if it is rank deficient
    # TODO: arima_constant should be considered to remove incompatible differences, rather than using the max
    xreg_c <- if(all(constant)) cbind(intercept = arima_constant(length(y), max(pdq$d), max(PDQ$D), period), xreg) else xreg
    qr <- qr(na.omit(xreg_c))
    num_regressors <- length(qr$qraux)
    if (qr$rank < num_regressors) {
      bad_regressors <- qr$pivot[(qr$rank + 1):num_regressors]
      # Ignore fixed coefficients for xreg
      bad_regressors <- bad_regressors[!(colnames(xreg_c)[bad_regressors] %in% union(names(fixed), names(specials$xreg[[1]]$fixed)))]
      # Offset for inclusion of intercept
      bad_regressors <- bad_regressors - all(constant)
      
      # Remove deficient regressors
      if(!is_empty(bad_regressors)){
        warn(sprintf(
          "Provided exogenous regressors are rank deficient, removing regressors: %s",
          paste("`", colnames(xreg)[bad_regressors], "`", sep = "", collapse = ", ")
        ))
        
        xreg <- xreg[, -bad_regressors, drop = FALSE]
      }
    }

    # Finally find residuals from regression in order
    # to estimate appropriate level of differencing
    j <- !is.na(x) & !is.na(rowSums(xreg))
    x[j] <- residuals(lm(x ~ xreg))
  }
  else {
    xreg <- NULL
  }

  diff <- function(x, differences, ...) {
    if (differences == 0) {
      return(x)
    }
    base::diff(x, differences = differences, ...)
  }

  # Choose seasonal differencing
  if (length(PDQ$D) > 1) {
    require_package("feasts")
    # Valid xregs

    if (!is.null(xreg)) {
      keep <- map_lgl(PDQ$D, function(.x) {
        if(.x*period >= nrow(xreg)) return(FALSE)
        diff_xreg <- diff(xreg, lag = period, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      PDQ$D <- PDQ$D[keep]
    }
    PDQ$D <- unname(feasts::unitroot_nsdiffs(stats::na.contiguous(x),
      alpha = unitroot_spec$nsdiffs_alpha,
      unitroot_fn = unitroot_spec$nsdiffs_pvalue,
      differences = PDQ$D, .period = period
    ))
  }
  x <- diff(x, lag = period, differences = PDQ$D)
  diff_xreg <- diff(xreg, lag = period, differences = PDQ$D)
  if (length(pdq$d) > 1) {
    require_package("feasts")

    # Valid xregs
    if (!is.null(xreg)) {
      keep <- map_lgl(pdq$d, function(.x) {
        if(.x >= nrow(diff_xreg)) return(FALSE)
        diff_xreg <- diff(diff_xreg, differences = .x)
        !any(apply(diff_xreg, 2, is.constant))
      })
      pdq$d <- pdq$d[keep]
    }

    pdq$d <- unname(feasts::unitroot_ndiffs(stats::na.contiguous(x),
      alpha = unitroot_spec$ndiffs_alpha,
      unitroot_fn = unitroot_spec$ndiffs_pvalue,
      differences = pdq$d
    ))
  }

  # Check number of differences selected
  if (length(PDQ$D) != 1) abort("Could not find appropriate number of seasonal differences.")
  if (length(pdq$d) != 1) abort("Could not find appropriate number of non-seasonal differences.")
  if (PDQ$D >= 2) {
    warn("Having more than one seasonal differences is not recommended. Please consider using only one seasonal difference.")
  } else if (PDQ$D + pdq$d > 2) {
    warn("Having 3 or more differencing operations is not recommended. Please consider reducing the total number of differences.")
  }

  # Find best model
  best <- NULL
  compare_arima <- function(p, d, q, P, D, Q, constant) {
    if (constant) {
      intercept <- arima_constant(length(y), d, D, period)
      xreg <- if (is.null(xreg)) {
        matrix(intercept, dimnames = list(NULL, "constant"))
      } else {
        cbind(xreg, intercept = intercept)
      }
    }
    
    if(is_empty(fixed)){
      fixed <- rep(NA_real_, p + q + P + Q + ncol(xreg)%||%0)
      fixed_nm <- c(sprintf("ar%i", seq_len(p)),sprintf("ma%i", seq_len(q)),sprintf("sar%i", seq_len(P)),sprintf("sma%i", seq_len(Q)),colnames(xreg))
      fixed_nm <- match(fixed_nm, names(special_fixed))
      fixed[which(!is.na(fixed_nm))] <- do.call("c", special_fixed[fixed_nm])
    }
    
    new <- wrap_arima(
      y,
      order = c(p, d, q),
      seasonal = list(order = c(P, D, Q), period = period),
      xreg = xreg, method = method, ..., include.mean = FALSE,
      fixed = fixed
    )

    if (!is.null(new)) {
      nstar <- length(y) - d - D * period
      npar <- length(new$coef[new$mask]) + 1
      if (method == "CSS") {
        new$aic <- offset + nstar * log(new$sigma2) + 2 * npar
      }
      if (!is.na(new$aic)) {
        new$bic <- new$aic + npar * (log(nstar) - 2)
        new$aicc <- new$aic + 2 * npar * (npar + 1) / (nstar - npar - 1)
      }
      else {
        new$aic <- new$bic <- new$aicc <- new$ic <- Inf
      }
      # Adjust residual variance to be unbiased
      new$sigma2 <- sum(new$residuals^2, na.rm = TRUE) / (nstar - npar + 1)

      # If automatically selecting a model
      if (NROW(model_opts) > 1) {
        # Check for unit roots
        minroot <- map_dbl(
          list(
            phi = -new$model$phi,
            theta = new$model$theta
          ),
          function(testvec) {
            k <- abs(testvec) > 1e-8
            if (any(k)) {
              last.nonzero <- max(which(k))
              testvec <- testvec[seq_len(last.nonzero)]
              min(abs(polyroot(c(1, testvec))))
            }
            else {
              2
            }
          }
        )

        if (isTRUE(min(minroot) < 1 + 1e-2)) { # Previously 1+1e-3
          new <- NULL
        } # Don't like this model
      }
    }
    if ((new[[ic]] %||% Inf) < (best[[ic]] %||% Inf)) {
      best <<- new
    }
    (new[[ic]] %||% Inf)
  }

  mostly_specified <- length(pdq$p) + length(pdq$d) + length(pdq$q) + length(PDQ$P) + length(PDQ$D) + length(PDQ$Q) == 6
  mostly_specified_msg <- "It looks like you're trying to fully specify your ARIMA model but have not said if a constant should be included.\nYou can include a constant using `ARIMA(y~1)` to the formula or exclude it by adding `ARIMA(y~0)`."
  model_opts <- expand.grid(p = pdq$p, d = pdq$d, q = pdq$q, P = PDQ$P, D = PDQ$D, Q = PDQ$Q, constant = constant)
  
  if (is.null(approximation)) {
    approximation <- ((length(x) > 150) || (period > 12)) && nrow(model_opts) > 1
  }
  if (approximation) {
    method <- "CSS"
    offset <- with(
      stats::arima(y,
                   order = c(0, pdq$d, 0), xreg = xreg,
                   include.mean = all(constant)
      ),
      -2 * loglik - NROW(data) * log(sigma2)
    )
  } else {
    method <- "CSS-ML"
  }
  
  if (NROW(model_opts) > 1) {
    model_opts <- filter(model_opts, !!enexpr(order_constraint))
    if (NROW(model_opts) == 0) {
      if (mostly_specified) warn(mostly_specified_msg)
      abort("There are no ARIMA models to choose from after imposing the `order_constraint`, please consider allowing more models.")
    }
    wrap_arima <- possibly(quietly(stats::arima), NULL)
  }
  else {
    wrap_arima <- stats::arima
  }

  if (any((model_opts$d + model_opts$D > 1) & model_opts$constant)) {
    warn("Model specification induces a quadratic or higher order polynomial trend. 
This is generally discouraged, consider removing the constant or reducing the number of differences.")
  }
  constant <- unique(model_opts$constant)

  if (stepwise) {
    # Prepare model comparison vector
    est_ic <- rep(NA_integer_, NROW(model_opts))
    best_ic <- Inf

    # Initial 4 models
    initial_opts <- list(
      start = c(pdq$p_init, pdq$d, pdq$q_init, PDQ$P_init, PDQ$D, PDQ$Q_init, constant[1]),
      null = c(0, pdq$d, 0, 0, PDQ$D, 0, constant[1]),
      ar = c(max(pdq$p) > 0, pdq$d, 0, max(PDQ$P) > 0, PDQ$D, 0, constant[1]),
      ma = c(0, pdq$d, max(pdq$q) > 0, 0, PDQ$D, max(PDQ$Q) > 0, constant[1])
    )
    step_order <- stats::na.omit(match(initial_opts, lapply(split(model_opts, seq_len(NROW(model_opts))), as.numeric)))
    initial <- TRUE

    # Stepwise search
    k <- 0
    while (NROW(model_opts[step_order, ]) > 0 && k < 94) {
      k <- k + 1

      # Evaluate model
      est_ic[step_order[1]] <- do.call(compare_arima, model_opts[step_order[1], ])

      if (greedy && !initial) {
        if (update_step <- est_ic[step_order[1]] < best_ic) {
          # Update best model and score
          best_ic <- est_ic[step_order[1]]
          current <- as.numeric(model_opts[step_order[1], ])
        }
      }
      else {
        if (update_step <- length(step_order) == 1) {
          best_ic <- min(est_ic, na.rm = TRUE)
          current <- as.numeric(model_opts[which.min(est_ic), ])
        }
      }

      if (update_step) {
        initial <- FALSE
        # Calculate new possible steps
        dist <- apply(model_opts, 1, function(x) sum((x - current)^2))
        step_order <- order(dist, model_opts$P, model_opts$Q, model_opts$p, model_opts$q)[seq_len(sum(dist <= 2))]
        step_order <- step_order[is.na(est_ic[step_order])]
      }
      else {
        # Move to next possible step
        step_order <- step_order[-1]
      }
    }
  }
  else {
    est_ic <- pmap_dbl(model_opts, compare_arima)
  }

  if (approximation && !is.null(best$arma)) {
    method <- "CSS-ML"
    best <- NULL
    step_order <- order(est_ic)[seq_len(sum(!is.na(est_ic)))]
    est_ic <- rep_len(Inf, length(est_ic)) # Ignore all approximate models until re-estimated
    for (mod_spec in step_order)
    {
      est_ic[mod_spec] <- do.call(compare_arima, model_opts[mod_spec, ])
      if (isTRUE(is.finite(est_ic[mod_spec]))) {
        break
      }
    }
  }

  if (is.null(best)) {
    if (mostly_specified) warn(mostly_specified_msg)
    abort("Could not find an appropriate ARIMA model.\nThis is likely because automatic selection does not select models with characteristic roots that may be numerically unstable.\nFor more details, refer to https://otexts.com/fpp3/arima-r.html#plotting-the-characteristic-roots")
  }

  # Compute ARMA roots
  ar <- best$model$phi
  ma <- best$model$theta
  arroot <- if (is_empty(ar) || !any(abs(ar) > 0)) cpl() else polyroot(c(1, -ar[seq_len(max(which(abs(ar) > 1e-8)))]))
  maroot <- if (is_empty(ma) || !any(abs(ma) > 0)) cpl() else polyroot(c(1, ma[seq_len(max(which(abs(ma) > 1e-8)))]))

  fit_coef <- coef(best)
  
  fit_se <- set_names(numeric(length(fit_coef)), names(fit_coef))
  fit_se[colnames(best$var.coef)] <- sqrt(diag(best$var.coef))
  
  if (is_empty(fit_se)) {
    fit_se <- NULL
  }
  else if (model_opts[which.min(est_ic), "constant"] && is.null(xreg)) {
    fit_coef["constant"] <- fit_coef["constant"] * (1 - sum(best$model$phi))
    fit_se["constant"] <- fit_se["constant"] * (1 - sum(best$model$phi))
  }

  # Compute regression residuals
  reg_resid <- as.numeric(y)
  if (model_opts[which.min(est_ic), "constant"]) {
    xreg <- cbind(xreg, constant = arima_constant(length(y), pdq$d, PDQ$D, period))
  }
  if (!is.null(xreg)) {
    reg_resid <- reg_resid - xreg %*% as.matrix(best$coef[(sum(best$arma[1:4]) + 1):length(best$coef)])
  }

  # Output model
  best_spec <- as_tibble(model_opts[which.min(est_ic),])
  best_spec[["period"]] <- period
  structure(
    list(
      par = tibble(
        term = names(fit_coef) %||% chr(), estimate = unname(fit_coef) %||% dbl(),
        std.error = fit_se %||% rep(NA, length(fit_coef))
      ) %>%
        mutate(
          statistic = !!sym("estimate") / !!sym("std.error"),
          p.value = 2 * stats::pt(abs(!!sym("statistic")),
            best$nobs,
            lower.tail = FALSE
          )
        ),
      est = tibble(
        .fitted = as.numeric(y - best$residuals),
        .resid = as.numeric(best$residuals),
        .regression_resid = reg_resid
      ),
      fit = tibble(
        sigma2 = best$sigma2,
        log_lik = best$loglik,
        AIC = best$aic, AICc = best$aicc, BIC = best$bic,
        ar_roots = list(arroot), ma_roots = list(maroot)
      ),
      spec = best_spec,
      model = best
    ),
    class = "ARIMA"
  )
}

specials_arima <- new_specials(
  pdq = function(p = 0:5, d = 0:2, q = 0:5,
                 p_init = 2, q_init = 2,
                 fixed = list()) {
    p <- p[p <= floor(NROW(self$data) / 3)]
    q <- q[q <= floor(NROW(self$data) / 3)]
    p_init <- p[which.min(abs(p - p_init))]
    q_init <- q[which.min(abs(q - q_init))]
    if(!all(grepl("^(ma|ar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for pdq() must begin with ar or ma, followed by a lag number.")
    }
    as.list(environment())
  },
  PDQ = function(P = 0:2, D = 0:1, Q = 0:2, period = NULL,
                 P_init = 1, Q_init = 1,
                 fixed = list()) {
    period <- get_frequencies(period, self$data, .auto = "smallest")
    if (period < 1) {
      abort("The seasonal period must be greater or equal to 1.")
    } else if (period == 1) {
      # Not seasonal
      P <- 0
      D <- 0
      Q <- 0
    }
    else {
      P <- P[P <= floor(NROW(self$data) / 3 / period)]
      Q <- Q[Q <= floor(NROW(self$data) / 3 / period)]
    }
    P_init <- P[which.min(abs(P - P_init))]
    Q_init <- Q[which.min(abs(Q - Q_init))]
    if(!all(grepl("^(sma|sar)\\d+", names(fixed)))){
      abort("The 'fixed' coefficients for PDQ() must begin with sar or sma, followed by a lag number.")
    }
    as.list(environment())
  },
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

    env$lag <- lag # Mask user defined lag to retain history when forecasting
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
  .required_specials = c("pdq", "PDQ"),
  .xreg_specials = names(common_xregs)
)

#' Estimate an ARIMA model
#'
#' Searches through the model space specified in the specials to identify the
#' best ARIMA model which has lowest AIC, AICc or BIC value. It is implemented
#' using [`stats::arima()`] and allows ARIMA models to be used in the fable
#' framework.
#'
#' @aliases report.ARIMA
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param stepwise Should stepwise be used?
#' @param greedy Should the stepwise search move to the next best option immediately?
#' @param approximation Should CSS (conditional sum of squares) be used during model selection? The default (`NULL`) will use the approximation if there are more than 150 observations or if the seasonal period is greater than 12.
#' @param order_constraint A logical predicate on the orders of `p`, `d`, `q`,
#' `P`, `D`, `Q` and `constant` to consider in the search. See "Specials" for 
#' the meaning of these terms.
#' @param unitroot_spec A specification of unit root tests to use in the
#' selection of `d` and `D`. See [`unitroot_options()`] for more details.
#' @param ... Further arguments for [`stats::arima()`]
#'
#' @section Parameterisation:
#'
#' The fable `ARIMA()` function uses an alternate parameterisation of constants
#' to [`stats::arima()`] and [`forecast::Arima()`]. While the parameterisations
#' are equivalent, the coefficients for the constant/mean will differ.
#'
#' In fable, the parameterisation used is:
#'
#' \deqn{(1-\phi_1B - \cdots - \phi_p B^p)(1-B)^d y_t = c + (1 + \theta_1 B + \cdots + \theta_q B^q)\varepsilon_t}
#'
#' In stats and forecast, an ARIMA model is parameterised as:
#'
#' \deqn{(1-\phi_1B - \cdots - \phi_p B^p)(y_t' - \mu) = (1 + \theta_1 B + \cdots + \theta_q B^q)\varepsilon_t}
#'
#' where \eqn{\mu} is the mean of \eqn{(1-B)^d y_t} and \eqn{c = \mu(1-\phi_1 - \cdots - \phi_p )}.
#'
#' @section Specials:
#' 
#' The _specials_ define the space over which `ARIMA` will search for the model that best fits the data. If the RHS of `formula` is left blank, the default search space is given by `pdq() + PDQ()`: that is, a model with candidate seasonal and nonseasonal terms, but no exogenous regressors. Note that a seasonal model requires at least 2 full seasons' worth of data; if this is not available, `ARIMA` will revert to a nonseasonal model with a warning.
#'
#' To specify a model fully (avoid automatic selection), the intercept and `pdq()/PDQ()` values must be specified: for example `formula = response ~ 1 + pdq(1, 1, 1) + PDQ(1, 0, 0)`.
#' 
#' \subsection{pdq}{
#' The `pdq` special is used to specify non-seasonal components of the model.
#' \preformatted{
#' pdq(p = 0:5, d = 0:2, q = 0:5,
#'     p_init = 2, q_init = 2, fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `p`      \tab The order of the non-seasonal auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `d`      \tab The order of integration for non-seasonal differencing. If multiple values are provided, one of the values will be selected via repeated KPSS tests. \cr
#'   `q`      \tab The order of the non-seasonal moving average (MA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `p_init` \tab If `stepwise = TRUE`, `p_init` provides the initial value for `p` for the stepwise search procedure. \cr
#'   `q_init` \tab If `stepwise = TRUE`, `q_init` provides the initial value for `q` for the stepwise search procedure. \cr
#'   `fixed`  \tab A named list of fixed parameters for coefficients. The names identify the coefficient, beginning with either `ar` or `ma`, and then followed by the lag order. For example, `fixed = list(ar1 = 0.3, ma2 = 0)`.
#' }
#' }
#'
#' \subsection{PDQ}{
#' The `PDQ` special is used to specify seasonal components of the model. To force a nonseasonal fit, specify `PDQ(0, 0, 0)` in the RHS of the model formula. Note that simply omitting `PDQ` from the formula will _not_ result in a nonseasonal fit.
#' \preformatted{
#' PDQ(P = 0:2, D = 0:1, Q = 0:2, period = NULL,
#'     P_init = 1, Q_init = 1, fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `P`      \tab The order of the seasonal auto-regressive (SAR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `D`      \tab The order of integration for seasonal differencing. If multiple values are provided, one of the values will be selected via repeated heuristic tests (based on strength of seasonality from an STL decomposition). \cr
#'   `Q`      \tab The order of the seasonal moving average (SMA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `period` \tab The periodic nature of the seasonality. This can be either a number indicating the number of observations in each seasonal period, or text to indicate the duration of the seasonal window (for example, annual seasonality would be "1 year").  \cr
#'   `P_init` \tab If `stepwise = TRUE`, `P_init` provides the initial value for `P` for the stepwise search procedure. \cr
#'   `Q_init` \tab If `stepwise = TRUE`, `Q_init` provides the initial value for `Q` for the stepwise search procedure. \cr
#'   `fixed`  \tab A named list of fixed parameters for coefficients. The names identify the coefficient, beginning with either `sar` or `sma`, and then followed by the lag order. For example, `fixed = list(sar1 = 0.1)`.
#' }
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an ARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`], where including `1` will add a constant and `0` or `-1` will remove the constant. If left out, the inclusion of a constant will be determined by minimising `ic`.
#'
#' \preformatted{
#' xreg(..., fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)\cr
#'   `fixed`    \tab A named list of fixed parameters for coefficients. The names identify the coefficient, and should match the name of the regressor. For example, `fixed = list(constant = 20)`.
#' }
#' }
#'
#' @seealso
#' [Forecasting: Principles and Practices, ARIMA models (chapter 9)](https://otexts.com/fpp3/arima.html)
#' [Forecasting: Principles and Practices, Dynamic regression models (chapter 10)](https://otexts.com/fpp3/dynamic.html)
#'
#' @return A model specification.
#'
#' @examples
#' # Manual ARIMA specification
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ 0 + pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   report()
#'
#' # Automatic ARIMA specification
#' library(tsibble)
#' library(dplyr)
#' tsibbledata::global_economy %>%
#'   filter(Country == "Australia") %>%
#'   model(ARIMA(log(GDP) ~ Population))
#' @importFrom stats model.matrix
#' @export
ARIMA <- function(formula, ic = c("aicc", "aic", "bic"), stepwise = TRUE, greedy = TRUE,
                  approximation = NULL,
                  order_constraint = p + q + P + Q <= 6 & (!constant | d + D < 2),
                  unitroot_spec = unitroot_options(), ...) {
  ic <- match.arg(ic)
  arima_model <- new_model_class("ARIMA",
    train = train_arima,
    specials = specials_arima, origin = NULL,
    check = all_tsbl_checks
  )
  new_model_definition(arima_model, !!enquo(formula),
    ic = ic,
    stepwise = stepwise, greedy = greedy,
    approximation = approximation,
    order_constraint = enexpr(order_constraint),
    unitroot_spec = unitroot_spec, ...
  )
}


#' Extract fitted values from a fable model
#'
#' Extracts the fitted values.
#'
#' @inheritParams forecast.ARIMA
#'
#' @return A vector of fitted values.
#'
#' @examples
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   fitted()
#' @export
fitted.ARIMA <- function(object, ...) {
  object$est[[".fitted"]]
}

#' Extract residuals values from a fable model
#'
#' Extracts the residuals.
#'
#' @inheritParams forecast.ARIMA
#' @param type The type of the residuals to extract.
#'
#' @return A vector of fitted residuals.
#'
#' @examples
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   residuals()
#' @export
residuals.ARIMA <- function(object, type = c("innovation", "regression"), ...) {
  type <- match.arg(type)
  if (type == "innovation") {
    object$est[[".resid"]]
  }
  else if (type == "regression") {
    object$est[[".regression_resid"]]
  }
  else {
    abort(sprintf('Residuals of `type = "%s"` are not supported by ARIMA models', type))
  }
}

#' Glance an ARIMA model
#'
#' Construct a single row summary of the ARIMA model.
#'
#' Contains the variance of residuals (`sigma2`), the log-likelihood (`log_lik`),
#' information criterion (`AIC`, `AICc`, `BIC`) and the characteristic roots
#' (`ar_roots` and `ma_roots`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   glance()
#' @export
glance.ARIMA <- function(x, ...) {
  x$fit
}

#' Tidy a fable model
#'
#' Returns the coefficients from the model in a `tibble` format.
#'
#' @inheritParams generics::tidy
#'
#' @return The model's coefficients in a `tibble`.
#'
#' @examples
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   tidy()
#' @export
tidy.ARIMA <- function(x, ...) {
  x$par
}

#' @export
report.ARIMA <- function(object, ...) {
  par <- rbind(tidy(object)$estimate, s.e. = tidy(object)$std.error)
  colnames(par) <- tidy(object)$term
  if (NCOL(par) > 0) {
    cat("\nCoefficients:\n")
    coef <- round(par, digits = 4)
    print.default(coef, print.gap = 2)
  }
  cat(
    "\nsigma^2 estimated as ", format(object$fit$sigma2, digits = 4),
    ":  log likelihood=", format(round(object$fit$log_lik, 2L)), "\n",
    sep = ""
  )

  cat("AIC=", format(round(object$fit$AIC, 2L)), sep = "")
  cat("   AICc=", format(round(object$fit$AICc, 2L)), sep = "")
  cat("   BIC=", format(round(object$fit$BIC, 2L)), "\n", sep = "")
}

#' Forecast a model from the fable package
#'
#' Produces forecasts from a trained model.
#'
#' @inheritParams fabletools::forecast
#' @param specials (passed by [`fabletools::forecast.mdl_df()`]).
#' @param bootstrap If `TRUE`, then forecast distributions are computed using simulation with resampled errors.
#' @param times The number of sample paths to use in estimating the forecast distribution when `boostrap = TRUE`.
#'
#' @importFrom stats formula residuals
#'
#' @return A list of forecasts.
#'
#' @examples
#' USAccDeaths %>%
#'   as_tsibble() %>%
#'   model(arima = ARIMA(log(value) ~ pdq(0, 1, 1) + PDQ(0, 1, 1))) %>%
#'   forecast()
#' @export
forecast.ARIMA <- function(object, new_data = NULL, specials = NULL,
                           bootstrap = FALSE, times = 5000, ...) {
  if (bootstrap) {
    sim <- map(seq_len(times), function(x) generate(object, new_data, specials, bootstrap = TRUE)[[".sim"]]) %>%
      transpose() %>%
      map(as.numeric)
    return(distributional::dist_sample(sim))
  }
  
  xreg <- specials$xreg[[1]]$xreg
  # Drop unused rank deficient xreg
  xreg <- xreg[,colnames(xreg) %in% names(object$model$coef), drop = FALSE]
  
  if (object$spec$constant) {
    intercept <- arima_constant(
      NROW(object$est) + NROW(new_data),
      object$spec$d, object$spec$D,
      object$spec$period
    )[NROW(object$est) + seq_len(NROW(new_data))]
    
    xreg <- if (is.null(xreg)) {
      matrix(intercept, dimnames = list(NULL, "constant"))
    } else {
      xreg <- cbind(xreg, intercept = intercept)
    }
  }
  # Produce predictions
  # Remove unnecessary warning for xreg
  object$model$call$xreg <- if (!is.null(xreg)) {
    expr(matrix(nrow = !!NROW(object$est), ncol = !!NCOL(xreg)))
  }
  else {
    NULL
  }
  fc <- predict(object$model, n.ahead = NROW(new_data), newxreg = xreg, ...)
  object$call$xreg <- NULL
  fc$pred <- as.numeric(fc$pred)
  fc$se <- as.numeric(fc$se)
  # Output forecasts
  distributional::dist_normal(fc$pred, fc$se)
}

#' @inherit generate.ETS
#' 
#' @examples 
#' fable_fit <- as_tsibble(USAccDeaths) %>%
#'   model(model = ARIMA(value ~ 0 + pdq(0,1,1) + PDQ(0,1,1)))
#' fable_fit %>% generate(times = 10)
#' 
#' @export
generate.ARIMA <- function(x, new_data, specials, bootstrap = FALSE, ...){
  # Get xreg
  h <- max(map_int(key_data(new_data)[[".rows"]], length))
  xreg <- specials$xreg[[1]]$xreg
  if(x$spec$constant){
    intercept <- arima_constant(NROW(x$est) + h,
                                x$spec$d, x$spec$D,
                                x$spec$period)[NROW(x$est) + seq_len(h)]
    
    xreg <- if(is.null(xreg)){
      matrix(intercept, dimnames = list(NULL, "constant"))
    } else {
      xreg <- cbind(xreg, intercept = intercept)
    }
  }
  
  if (!is.null(xreg)) {
    narma <- sum(x$model$arma[1L:4L])
    coef <- x$model$coef
    coef <- coef[(narma + 1):length(coef)]
    xm <- drop(xreg %*% coef)
  } else {
    xm <- rep(0, nrow(new_data))
  }
  
  if(!(".innov" %in% names(new_data))){
    if(bootstrap){
      res <- residuals(x)
      new_data$.innov <- sample(na.omit(res) - mean(res, na.rm = TRUE),
                                nrow(new_data), replace = TRUE)
    }
    else{
      new_data$.innov <- stats::rnorm(nrow(new_data), sd = sqrt(x$fit$sigma2))
    }
  }
  
  new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := conditional_arima_sim(x$model, x$est$.regression_resid, !!sym(".innov"))) %>% 
    dplyr::ungroup() %>% 
    mutate(".sim" := as.numeric(!!sym(".sim") + xm))
}

# Version of stats::arima.sim which conditions on past observations
#' @importFrom stats diffinv
conditional_arima_sim <- function(object, x, e){
  theta <- object$model$theta
  phi <- object$model$phi
  m <- object$arma[5L]
  d <- object$arma[6L]
  D <- object$arma[7L]
  
  past <- residuals(object)
  y <- ts(c(past, e), frequency = object$arma[5])
  # MA filtering
  if (any(theta != 0)) {
    y <- stats::filter(y, c(1, theta), method = "convolution", sides = 1L)
    y[seq_along(theta)] <- 0
  }
  y <- y[seq(length(past) + 1, length(y))]
  
  # AR filtering
  if (any(phi != 0)) {
    dx <- x
    if (d > 0){
      dx <- diff(dx, lag = 1, differences = d)
    }
    if (D > 0){
      dx <- diff(dx, lag = m, differences = D)
    }
    
    y <- stats::filter(y, phi, method = "recursive", init = rev(dx)[1:length(phi)])
    # TODO: fill missing init data with zero
  }
  
  # Undo differences
  if (d > 0){
    # Regular undifferencing
    if (D == 0) {
      # if there is no seasonal differencing
      y <- diffinv(y, differences = d, xi = x[length(x) - (d:1) + 1])[-(1:d)]
    } else {
      # if there is seasonal differencing
      dx <- diff(x, lag = m, differences = D)
      dxi <- dx[(length(dx) - D):length(dx)]
      y <- diffinv(y, differences = d, xi = dxi[length(dxi) - (d:1) + 1])[-(1:d)]
    }
  }
  if (D > 0) {
    # Seasonal undifferencing
    i <- length(x) - D * m + 1
    xi <- x[i:length(x)]
    y <- diffinv(y, lag = m, differences = D, xi = xi)[-(1:length(xi))]
  }
  return(y)
}

#' Refit an ARIMA model
#'
#' Applies a fitted ARIMA model to a new dataset.
#'
#' @inheritParams forecast.ARIMA
#' @param reestimate If `TRUE`, the coefficients for the fitted model will be re-estimated to suit the new data.
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(ARIMA(value ~ 1 + pdq(2, 0, 0) + PDQ(2, 1, 0)))
#'
#' report(fit)
#'
#' fit %>%
#'   refit(lung_deaths_female) %>%
#'   report()
#' @return A refitted model.
#'
#' @importFrom stats formula residuals
#' @export
refit.ARIMA <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  # Update data for re-evaluation
  specials$pdq[[1]][c("p", "d", "q", "p_init", "q_init")] <- 
    as.list(object$spec[c("p", "d", "q", "p", "q")])
  specials$PDQ[[1]][c("P", "D", "Q", "period", "P_init", "Q_init")] <-
    as.list(object$spec[c("P", "D", "Q", "period", "P", "Q")])
  
  if (reestimate) {
    return(train_arima(new_data, specials, ...))
  }

  out <- train_arima(new_data, specials, fixed = object$model$coef, ...)
  out$par <- object$par
  out
}

#' Interpolate missing values from a fable model
#'
#' Applies a model specific estimation technique to predict the values of missing values in a `tsibble`, and replace them.
#'
#' @inheritParams forecast.ARIMA
#'
#' @return A tibble of the same dimension of `new_data` with missing values interpolated.
#'
#' @examples
#' library(tsibbledata)
#'
#' olympic_running %>%
#'   model(arima = ARIMA(Time ~ trend())) %>%
#'   interpolate(olympic_running)
#' @importFrom stats formula residuals
#' @export
interpolate.ARIMA <- function(object, new_data, specials, ...) {
  # Get missing values
  y <- unclass(new_data)[[measured_vars(new_data)]]
  miss_val <- which(is.na(y))
  object <- refit(object, new_data, specials, ...)$model$model
  fits <- stats::KalmanSmooth(y, object)$smooth[miss_val, , drop = FALSE] %*% as.matrix(object$Z)

  # Update data
  i <- (miss_val - 1) %% NROW(new_data) + 1
  j <- (miss_val - 1) %/% NROW(new_data) + 1
  idx_pos <- match(index_var(new_data), colnames(new_data))
  j <- ifelse(j >= idx_pos, j + 1, j)
  pos <- split(i, j)
  for (i in seq_along(pos)) {
    new_data[[as.numeric(names(pos)[i])]][pos[[i]]] <- fits
  }
  new_data
  new_data
}

#' @export
model_sum.ARIMA <- function(x) {
  out <- sprintf(
    "ARIMA(%i,%i,%i)%s",
    x$spec$p, x$spec$d, x$spec$q,
    if (any(x$spec[c("P", "D", "Q")] > 0)) {
      sprintf(
        "(%i,%i,%i)[%i]",
        x$spec$P, x$spec$D, x$spec$Q, x$spec$period
      )
    } else {
      ""
    }
  )
  if (NROW(x$par) > sum(x$spec[c("p", "q", "P", "Q", "constant")])) {
    out <- sprintf("LM w/ %s errors", out)
  }
  else if (x$spec$constant) {
    out <- sprintf(
      "%s w/ %s", out,
      if (x$spec$d + x$spec$D == 0) {
        "mean"
      } else if (x$spec$d + x$spec$D == 1) {
        "drift"
      } else {
        "poly"
      }
    )
  }
  out
}

arima_constant <- function(n, d, D, period) {
  constant <- rep(1, n)
  if (d > 0) {
    constant <- stats::diffinv(constant, differences = d, xi = rep(1, d))[seq_len(n)]
  }
  if (D > 0) {
    constant <- stats::diffinv(constant, lag = period, differences = D, xi = rep(1, period * D))[seq_len(n)]
  }
  constant
}

#' Options for the unit root tests for order of integration
#' 
#' By default, a kpss test (via [`feasts::unitroot_kpss()`]) will be performed 
#' for testing the required first order differences, and a test of the seasonal 
#' strength (via [`feasts::feat_stl()`] seasonal_strength) being above the 0.64 
#' threshold is used for determining seasonal required differences.
#'
#' @param ndiffs_alpha,nsdiffs_alpha The level for the test specified in the `pval` functions As long as `pval < alpha`, differences will be added.
#' @param ndiffs_pvalue,nsdiffs_pvalue A function (or lambda expression) which returns the probability of the . As long as `pval < alpha`, differences will be added.
#'
#' For the function for the seasonal p-value, the seasonal period will be provided as the `.period` argument to this function.
#' A vector of data to test is available as `.` or `.x`.
#'
#' @return A list of parameters
#'
#' @export
unitroot_options <- function(ndiffs_alpha = 0.05, nsdiffs_alpha = 0.05,
                             ndiffs_pvalue = ~ feasts::unitroot_kpss(.)["kpss_pvalue"],
                             nsdiffs_pvalue = ur_seasonal_strength(0.64)) {
  as.list(environment())
}

globalVariables(".period")

ur_seasonal_strength <- function(threshold = 0.64){
  function(x){
    features <- feasts::feat_stl(x, .period)
    seas_strength <- grepl("^seasonal_strength_", names(features))
    if(!any(seas_strength)){
      FALSE
    } else{
      features[seas_strength] < threshold
    }
  }
}