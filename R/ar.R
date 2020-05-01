#' Estimate a AR model
#'
#' Searches through the vector of lag orders to find the best AR model which
#' has lowest AIC, AICc or BIC value. It is implemented using OLS, and behaves
#' comparably to [`stats::ar.ols()`].
#'
#' Exogenous regressors and [`common_xregs`] can be specified in the model
#' formula.
#'
#' @aliases report.AR
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param ... Further arguments for arima
#'
#' @section Specials:
#'
#' \subsection{pdq}{
#' The `order` special is used to specify the lag order for the auto-regression.
#' \preformatted{
#' order(p = 0:15, fixed = list())
#' }
#'
#' \tabular{ll}{
#'   `p`        \tab The order of the auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen.\cr
#'   `fixed`    \tab A named list of fixed parameters for coefficients. The names identify the coefficient, beginning with `ar`, and then followed by the lag order. For example, `fixed = list(ar1 = 0.3, ar3 = 0)`.
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
#' @return A model specification.
#'
#' @seealso
#' [Forecasting: Principles and Practices, Vector autoregressions (section 11.2)](https://otexts.com/fpp3/AR.html)
#'
#' @examples
#' luteinizing_hormones <- as_tsibble(lh)
#' fit <- luteinizing_hormones %>%
#'   model(AR(value ~ order(3)))
#'
#' report(fit)
#'
#' fit %>%
#'   forecast() %>%
#'   autoplot(luteinizing_hormones)
#' @export
AR <- function(formula, ic = c("aicc", "aic", "bic"), ...) {
  ic <- match.arg(ic)
  ar_model <- new_model_class("AR",
                                 train = train_ar,
                                 specials = specials_ar,
                                 origin = NULL,
                                 check = all_tsbl_checks
  )
  new_model_definition(ar_model, !!enquo(formula), ic = ic, ...)
}

specials_ar <- new_specials(
  order = function(p = 0:15, fixed = list()) {
    if (any(p < 0)) {
      warn("The AR order must be non-negative. Only non-negative orders will be considered.")
      p <- p[p >= 0]
    }
    list(p = p, fixed = fixed)
  },
  common_xregs,
  xreg = function(..., fixed = list()) {
    dots <- enexprs(...)
    env <- map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
    env <- if (!is_empty(env)) get_env(env[[1]]) else base_env()
    
    constants <- map_lgl(dots, inherits, "numeric")
    constant_given <- any(map_lgl(dots[constants], `%in%`, -1:1))
    
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    list(
      constant = if (constant_given) as.logical(terms(xreg) %@% "intercept") else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg,
      fixed = fixed
    )
  },
  .required_specials = c("order", "xreg"),
  .xreg_specials = names(common_xregs)
)

#' @importFrom stats ts
train_ar <- function(.data, specials, ic, ...) {
  # Get args
  p <- specials$order[[1]]$p
  
  # Get response variables
  y <- unclass(.data)[[measured_vars(.data)]]
  
  # Get xreg
  constant <- specials$xreg[[1]]$constant %||% c(TRUE, FALSE)
  xreg <- specials$xreg[[1]]$xreg
  
  fixed <- c(specials$order[[1]]$fixed, specials$xreg[[1]]$fixed)
  
  # Choose best model
  reduce(transpose(expand.grid(p = p, constant = constant)),
         function(best, args) {
           new <- estimate_ar(y, args$p, xreg, args$constant, fixed)
           if ((new[[ic]] %||% Inf) < (best[[ic]] %||% Inf)) {
             best <- new
           }
           best
         },
         .init = NULL
  )
}

# Adapted and generalised from stats::ar.ols
estimate_ar <- function(x, p, xreg, constant, fixed) {
  if (is.null(xreg)) {
    xreg <- matrix(nrow = length(x), ncol = 0)
  }
  if (constant) {
    xreg <- cbind(constant = rep.int(1, length(x)), xreg)
  }
  
  # scale
  x_sd <- sd(x)
  x <- x/x_sd
    
  par <- c(colnames(xreg), sprintf("ar%i", seq_len(p)))
  coef <- set_names(map_dbl(fixed[par], `%||%`, NA_real_), par)
  
  y <- stats::embed(x, p + 1L)
  X <- cbind(xreg[(p+1):nrow(xreg),,drop=FALSE], y[,-1,drop=FALSE])
  Y <- y[,1]
  Y_est <- t(Y - X[,!is.na(coef),drop=FALSE]%*%coef[!is.na(coef)])
  X_est <- X[,is.na(coef),drop=FALSE]
  
  nobs <- length(x)
  npar <- ncol(X_est)
  nr <- nrow(X_est)
  
  XX <- t(X_est) %*% X_est
  rank <- qr(XX)$rank
  if (rank != nrow(XX)) {
    warning(paste("model order: ", p, "singularities in the computation of the projection matrix", 
                  "results are only valid up to model order", 
                  p - 1L), domain = NA)
    return(NULL)
  }
  P <- if (ncol(XX) > 0) 
    solve(XX)
  else XX
  coef[coef_est <- is.na(coef)] <- Y_est %*% X_est %*% P
  YH <- drop(coef %*% t(X))
  E <- Y - YH
  varE <- tcrossprod(t(E))/nr
  varA <- kronecker(P, varE)
  
  coef_se <- numeric(length(par))
  coef_se[coef_est] <- if (ncol(varA) > 0) sqrt(diag(varA)) else numeric()
  
  aic <- nobs * log(det(varE)) + 2 * npar
  bic <- aic + npar * (log(nobs) - 2)
  aicc <- aic + 2 * npar * (npar + 1) / (nobs - npar - 1)
  
  # rescale
  coef[seq_len(ncol(xreg))] <- coef[seq_len(ncol(xreg))]*x_sd
  coef_se[seq_len(ncol(xreg))] <- coef_se[seq_len(ncol(xreg))]*x_sd
  varE <- varE * x_sd^2
  E <- E * x_sd
  x <- x * x_sd
  
  # Output model
  structure(
    list(
      coef = coef,
      coef.se = coef_se,
      fits = c(rep.int(NA_real_, p), YH),
      resid = c(rep.int(NA_real_, p), E),
      last = x[(length(E)+1):length(x)],
      sigma2 = drop(varE),
      aic = aic,
      bic = bic,
      aicc = aicc,
      p = p, 
      constant = constant
    ),
    class = "AR"
  )
}

#' @inherit forecast.ARIMA
#' @examples
#' as_tsibble(lh) %>%
#'   model(AR(value ~ order(3))) %>%
#'   forecast()
#' @export
forecast.AR <- function(object, new_data = NULL, specials = NULL,
                         bootstrap = FALSE, times = 5000, ...) {
  if (bootstrap) {
    abort("Bootstrapped forecasts for ARs are not yet implemented.")
  }
  
  h <- NROW(new_data)
  x <- c(object$last, rep(NA_real_, h))
  n <- length(object$last)
  p <- object$p
  coef <- object$coef
  # Get xreg
  xreg <- specials$xreg[[1]]$xreg
  if (object$constant) {
    xreg <- cbind(constant = rep(1, h), xreg)
  }
  
  # Predict
  nx <- length(coef) - p
  if (!is.null(xreg)) {
    xcoef <- coef[seq_len(nx)]
    xm <- drop(xreg %*% xcoef)
  } else {
    xm <- rep(0, nrow(new_data))
  }
  
  ar <- coef[nx + seq_len(p)]
  for (i in seq_len(h)) x[p + i] <- sum(ar * x[n + i - seq_len(p)]) + xm[i]
  
  fc <- x[p + seq_len(h)]
  psi <- ar_se(ar, h)
  se <- sqrt(object$sigma2 * cumsum(c(1, psi^2)))
  
  # Output forecasts
  distributional::dist_normal(fc, se)
}

#' Refit an AR model
#'
#' Applies a fitted AR model to a new dataset.
#'
#' @inheritParams forecast.AR
#' @param reestimate If `TRUE`, the coefficients for the fitted model will be re-estimated to suit the new data.
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(AR(value ~ 1 + order(10)))
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
refit.AR <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  y <- unclass(new_data)[[measured_vars(new_data)]]
  
  fixed <- if (reestimate) {
    c(specials$order[[1]]$fixed, specials$xreg[[1]]$fixed)
  } else {
    as.list(object$coef)
  }
  
  estimate_ar(y, object$p, specials$xreg[[1]]$xreg, object$constant, fixed)
}

ar_se <- function(phi, h){
  p <- length(drop(phi))
  npsi <- h + p
  psi <- numeric(npsi)
  psi[seq_len(p)] <- phi[seq_len(p)]
  for(i in seq_len(h))
    for(j in seq_len(p)) psi[i + j] = psi[i + j] + phi[j] * psi[i]
  psi[seq_len(h-1)]
  # vars <- cumsum(c(1, psi^2))
  # sqrt(object$var.pred * vars)[seq_len(h)]
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' as_tsibble(lh) %>%
#'   model(AR(value ~ order(3))) %>%
#'   fitted()
#' @export
fitted.AR <- function(object, ...) {
  object$fits
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' as_tsibble(lh) %>%
#'   model(AR(value ~ order(3))) %>%
#'   residuals()
#' @export
residuals.AR <- function(object, type = c("innovation", "regression"), ...) {
  type <- match.arg(type)
  if (type == "innovation") {
    object$resid
  }
  else if (type == "regression") {
    object$reg_resid
  }
}

#' @export
model_sum.AR <- function(x) {
  sprintf("AR(%s)%s", x$p, ifelse(x$constant, " w/ mean", ""))
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' as_tsibble(lh) %>%
#'   model(AR(value ~ order(3))) %>%
#'   tidy()
#' @export
tidy.AR <- function(x) {
  out <- tibble::enframe(drop(x$coef), "term", "estimate")
  out$std.error <- x$coef.se 
  out$statistic <- out$estimate / out$std.error
  out$p.value <- 2*stats::pt(abs(out$statistic), length(x$resid) - length(x$coef), lower.tail = FALSE)
  out
}

#' Glance a AR
#'
#' Construct a single row summary of the AR model.
#'
#' Contains the variance of residuals (`sigma2`), the log-likelihood (`log_lik`),
#' and information criterion (`AIC`, `AICc`, `BIC`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' as_tsibble(lh) %>%
#'   model(AR(value ~ order(3))) %>%
#'   glance()
#' @export
glance.AR <- function(x, ...) {
  tibble(sigma2 = x$sigma2, AIC = x$aic, AICc = x$aicc, BIC = x$bic, dof = length(x$resid) - length(x$coef))
}

#' @export
report.AR <- function(object, ...) {
  par <- rbind(tidy(object)$estimate)
  colnames(par) <- tidy(object)$term
  rownames(par) <- ""
  if (NCOL(par) > 0) {
    cat("\nCoefficients:\n")
    coef <- round(par, digits = 4)
    print.default(coef, print.gap = 2)
  }
  cat(
    "\nsigma^2 estimated as ", format(object$sigma2, digits = 4),
    # ":  log likelihood=", format(round(object$fit$log_lik, 2L)),
    "\n",
    sep = ""
  )
  cat(
    sprintf(
      "AIC = %s\tAICc = %s\tBIC = %s",
      format(round(object$aic, 2L)),
      format(round(object$aicc, 2L)),
      format(round(object$bic, 2L))
    )
  )
}
