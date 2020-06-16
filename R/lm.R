lm_glance_measures <- function(fit) {
  # Set up fit measures
  res <- fit$residuals[complete.cases(fit$residuals), , drop = FALSE]
  fits <- fit$fitted[complete.cases(fit$fitted), , drop = FALSE]
  n <- length(res)
  rdf <- fit$df.residual
  edf <- n - rdf
  rank <- fit$rank

  coef <- fit$coefficients
  intercept <- "(Intercept)" %in% rownames(coef)

  mss <- sum((fits - intercept * mean(fits))^2)
  rss <- sum(res^2)
  resvar <- rss / rdf

  if (NROW(coef) - intercept == 0) {
    r.squared <- adj.r.squared <- 0
    fstatistic <- NA
    p.value <- NA
  }
  else {
    r.squared <- mss / (mss + rss)
    adj.r.squared <- 1 - (1 - r.squared) * ((n - intercept) / rdf)
    fstatistic <- (mss / (rank - intercept)) / resvar
    p.value <- stats::pf(fstatistic, rank - intercept, rdf, lower.tail = FALSE)
  }

  influence <- stats::lm.influence(fit)

  dev <- n * log(rss / n)
  aic <- dev + 2 * edf + 2
  k <- edf - 1

  loglik <- 0.5 * (-n * (log(2 * pi) + 1 - log(n) + log(rss)))

  list(
    r_squared = r.squared, adj_r_squared = adj.r.squared,
    sigma2 = resvar, statistic = fstatistic,
    p_value = p.value, df = edf, log_lik = loglik,
    AIC = aic, AICc = aic + 2 * (k + 2) * (k + 3) / (n - k - 3),
    BIC = aic + (k + 2) * (log(n) - 2),
    CV = mean((res / (1 - influence$hat))^2, na.rm = TRUE),
    deviance = rss, df.residual = rdf, rank = rank
  )
}

lm_tidy_measures <- function(fit) {
  tibble(
    term = names(coef(fit))[!is.na(coef(fit))] %||% chr(),
    !!!as_tibble(`colnames<-`(coef(summary(fit)), c("estimate", "std.error", "statistic", "p.value")))
  )
}

train_tslm <- function(.data, specials, ...) {
  y <- invoke(cbind, unclass(.data)[measured_vars(.data)])
  xreg <- specials$xreg[[1]]

  keep <- complete.cases(xreg) & complete.cases(y)
  fit <- stats::lm.fit(xreg[keep, , drop = FALSE], y[keep, , drop = FALSE])
  resid <- matrix(nrow = nrow(y), ncol = ncol(y))
  resid[keep, ] <- as.matrix(fit$residuals)
  fit$residuals <- resid
  fit$fitted.values <- y - fit$residuals

  if (is_empty(fit$coefficients)) {
    fit$coefficients <- matrix(nrow = 0, ncol = NCOL(y))
  }
  else {
    fit$coefficients <- as.matrix(fit$coefficients)
  }
  colnames(fit$coefficients) <- colnames(y)
  
  # Remove unused structure
  fit$effects <- NULL
  fit$sigma2 <- sum(resid^2, na.rm = TRUE)/fit$df.residual 
  
  structure(fit, class = "TSLM")
}

specials_tslm <- new_specials(
  common_xregs,
  xreg = function(...) {
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(enexprs(...), function(.x, .y) call2("+", .x, .y))
    )
    env <- parent.frame()
    if (!exists("list", env)) env <- base_env()

    env$lag <- lag # Mask user defined lag to retain history when forecasting
    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    mm <- model.matrix(terms(xreg), xreg)
    if (NROW(mm) == 0 && identical(colnames(mm), "(Intercept)")) {
      return(matrix(data = 1, nrow = NROW(self$data), dimnames = list(NULL, "(Intercept)")))
    }
    mm
  },
  .required_specials = "xreg",
  .xreg_specials = names(common_xregs),
)

#' Fit a linear model with time series components
#'
#' The model formula will be handled using [`stats::model.matrix()`], and so
#' the the same approach to include interactions in [`stats::lm()`] applies when
#' specifying the `formula`. In addition to [`stats::lm()`], it is possible to
#' include [`common_xregs`] in the model formula, such as `trend()`, `season()`,
#' and `fourier()`.
#'
#' @aliases report.TSLM
#'
#' @param formula Model specification.
#'
#' @section Specials:
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an ARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#' \preformatted{
#' xreg(...)
#' }
#'
#' \tabular{ll}{
#'   `...`      \tab Bare expressions for the exogenous regressors (such as `log(x)`)
#' }
#' }
#'
#' @return A model specification.
#'
#' @seealso
#' [`stats::lm()`], [`stats::model.matrix()`]
#' [Forecasting: Principles and Practices, Time series regression models (chapter 6)](https://otexts.com/fpp3/regression.html)
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season()))
#'
#' library(tsibbledata)
#' olympic_running %>%
#'   model(TSLM(Time ~ trend())) %>%
#'   interpolate(olympic_running)
#' @export
TSLM <- function(formula) {
  tslm_model <- new_model_class("TSLM",
    train = train_tslm,
    specials = specials_tslm, origin = NULL
  )
  new_model_definition(tslm_model, !!enquo(formula))
}

#' @inherit fitted.ARIMA
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   fitted()
#' @export
fitted.TSLM <- function(object, ...) {
  object$fitted
}

#' @inherit residuals.ARIMA
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   residuals()
#' @export
residuals.TSLM <- function(object, ...) {
  object$residuals
}

#' Glance a TSLM
#'
#' Construct a single row summary of the TSLM model.
#'
#' Contains the R squared (`r_squared`), variance of residuals (`sigma2`),
#' the log-likelihood (`log_lik`), and information criterion (`AIC`, `AICc`, `BIC`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   glance()
#' @export
glance.TSLM <- function(x, ...) {
  as_tibble(lm_glance_measures(x))
}

#' @inherit tidy.ARIMA
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   tidy()
#' @export
tidy.TSLM <- function(x, ...) {
  rdf <- x$df.residual
  coef <- x$coefficients
  rank <- x$rank

  R <- chol2inv(x$qr$qr[seq_len(rank), seq_len(rank), drop = FALSE])

  se <- rep(NA_real_, nrow(coef))
  se[!is.na(coef)] <- sqrt(diag(R) * x$sigma2) # map(resvar, function(resvar) sqrt(diag(R) * resvar)) #for multiple response tslm

  out <- dplyr::as_tibble(x$coefficients, rownames = "term") %>%
    tidyr::gather(".response", "estimate", !!!syms(colnames(coef)))
  if (NCOL(coef) == 1) out[[".response"]] <- NULL
  out %>%
    mutate(
      std.error = unlist(se),
      statistic = !!sym("estimate") / !!sym("std.error"),
      p.value = 2 * stats::pt(abs(!!sym("statistic")), rdf, lower.tail = FALSE)
    )
}

#' @export
report.TSLM <- function(object, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nResiduals:\n")
  glance <- glance(object)
  intercept <- "(Intercept)" %in% rownames(object$coef)

  rdf <- glance$df.residual
  if (rdf > 5L) {
    res <- residuals(object)
    res_qt <- zapsmall(stats::quantile(res, na.rm = TRUE))
    names(res_qt) <- c("Min", "1Q", "Median", "3Q", "Max")
    print(res_qt, digits = digits, ...)
  }

  cat("\nCoefficients:\n")
  coef <- tidy(object)
  coef_mat <- as.matrix(coef[ncol(coef) - c(3:0)])
  colnames(coef_mat) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(coef_mat) <- coef$term
  stats::printCoefmat(coef_mat,
    digits = digits,
    signif.stars = getOption("show.signif.stars"), ...
  )

  cat(sprintf(
    "\nResidual standard error: %s on %s degrees of freedom\n",
    format(signif(sqrt(glance$sigma2), digits)), rdf
  ))
  if (!is.na(glance$statistic)) {
    cat(sprintf(
      "Multiple R-squared: %s,\tAdjusted R-squared: %s\nF-statistic: %s on %s and %s DF, p-value: %s\n",
      formatC(glance$r_squared, digits = digits),
      formatC(glance$adj_r_squared, digits = digits),
      formatC(glance$statistic, digits = digits),
      format(glance$rank - intercept), format(rdf),
      format.pval(glance$p_value)
    ))
  }
  invisible(object)
}

#' @inherit forecast.ARIMA
#' @importFrom stats predict
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   forecast()
#' @export
forecast.TSLM <- function(object, new_data, specials = NULL, bootstrap = FALSE,
                          times = 5000, ...) {
  coef <- object$coefficients
  rank <- object$rank
  qr <- object$qr
  piv <- qr$pivot[seq_len(rank)]

  # Get xreg
  xreg <- specials$xreg[[1]]

  if (rank < ncol(xreg)) {
    warn("prediction from a rank-deficient fit may be misleading")
  }

  # Intervals
  if (bootstrap) { # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x) {
      generate(object, new_data, specials, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose() %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    distributional::dist_sample(sim)
  } else {
    fc <- xreg[, piv, drop = FALSE] %*% coef[piv]
    resvar <- object$sigma2

    if (rank > 0) {
      XRinv <- xreg[, piv] %*% qr.solve(qr.R(qr)[seq_len(rank), seq_len(rank)])
      ip <- drop(XRinv^2 %*% rep(resvar, rank))
    }
    else {
      ip <- rep(0, length(fc))
    }

    se <- sqrt(ip + resvar)
    distributional::dist_normal(drop(fc), drop(se))
  }
}

#' @inherit generate.ETS
#'
#' @examples
#' as_tsibble(USAccDeaths) %>%
#'   model(lm = TSLM(log(value) ~ trend() + season())) %>%
#'   generate()
#' @export
generate.TSLM <- function(x, new_data, specials, bootstrap = FALSE, ...) {
  # Get xreg
  xreg <- specials$xreg[[1]]

  coef <- x$coefficients
  piv <- x$qr$pivot[seq_len(x$rank)]
  pred <- xreg[, piv, drop = FALSE] %*% coef[piv]

  if (!(".innov" %in% names(new_data))) {
    if (bootstrap) {
      res <- residuals(x)
      new_data$.innov <- sample(na.omit(res) - mean(res, na.rm = TRUE),
        NROW(new_data),
        replace = TRUE
      )
    }
    else {
      vars <- x$sigma2 / x$df.residual
      new_data$.innov <- stats::rnorm(length(pred), sd = sqrt(vars))
    }
  }

  transmute(new_data, .sim = pred + !!sym(".innov"))
}

#' @inherit interpolate.ARIMA
#'
#' @examples
#' library(tsibbledata)
#'
#' olympic_running %>%
#'   model(lm = TSLM(Time ~ trend())) %>%
#'   interpolate(olympic_running)
#' @export
interpolate.TSLM <- function(object, new_data, specials, ...) {
  # Get inputs
  miss_val <- which(is.na(new_data[measured_vars(new_data)]))
  xreg <- specials$xreg[[1]]

  # Make predictions
  coef <- object$coefficients
  piv <- object$qr$pivot[seq_len(object$rank)]
  pred <- xreg[miss_val, piv, drop = FALSE] %*% coef[piv]

  # Update data
  i <- (miss_val - 1) %% NROW(new_data) + 1
  j <- (miss_val - 1) %/% NROW(new_data) + 1
  idx_pos <- match(index_var(new_data), colnames(new_data))
  j <- ifelse(j >= idx_pos, j + 1, j)
  pos <- split(i, j)
  for (i in seq_along(pos)) {
    new_data[[as.numeric(names(pos)[i])]][pos[[i]]] <- pred
  }
  new_data
}

#' Refit a `TSLM`
#'
#' Applies a fitted `TSLM` to a new dataset.
#'
#' @inheritParams refit.ARIMA
#'
#' @examples
#' lung_deaths_male <- as_tsibble(mdeaths)
#' lung_deaths_female <- as_tsibble(fdeaths)
#'
#' fit <- lung_deaths_male %>%
#'   model(TSLM(value ~ trend() + season()))
#'
#' report(fit)
#'
#' fit %>%
#'   refit(lung_deaths_female) %>%
#'   report()
#' @export
refit.TSLM <- function(object, new_data, specials = NULL, reestimate = FALSE, ...) {
  # Update data for re-evaluation
  if (reestimate) {
    return(train_tslm(new_data, specials, ...))
  }

  # Get inputs
  y <- invoke(cbind, unclass(new_data)[measured_vars(new_data)])
  xreg <- specials$xreg[[1]]

  fit <- object
  coef <- object$coefficients
  fit$qr <- qr(xreg)
  piv <- fit$qr$pivot[seq_len(fit$rank)]
  fit$fitted.values <- xreg[, piv, drop = FALSE] %*% coef[piv]
  fit$residuals <- y - fit$fitted.values
  fit$sigma2 <- sum(fit$residuals^2, na.rm = TRUE)/fit$df.residual 

  structure(fit, class = "TSLM")
}

#' @export
model_sum.TSLM <- function(x) {
  "TSLM"
}
