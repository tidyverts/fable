train_varima <- function(.data, specials, ic, ...) {
  # Get response variables
  y <- invoke(cbind, unclass(.data)[measured_vars(.data)])
  
  if(any(colnames(specials$xreg[[1]]) != "(Intercept)")) {
    stop("Exogenous regressors for VARIMA are not yet supported.")
  }
  
  p <- specials$pdq[[1]]$p
  d <- specials$pdq[[1]]$d
  q <- specials$pdq[[1]]$q
  
  require_package("MTS")
  utils::capture.output(
    fit <- MTS::VARMA(
      if(d > 0) diff(y, differences = d) else y,
      p = p, q = q,
      include.mean = "(Intercept)" %in% colnames(specials$xreg[[1]])
    )
  )
  dimnames(fit$Sigma) <- list(colnames(y), colnames(y))
  fit$y_start <- y[seq_len(d),,drop = FALSE]
  fit$y_end <- y[NROW(y) - d + seq_len(d),,drop = FALSE]
  structure(fit, class="VARIMA")
}

specials_varima <- new_specials(
  pdq = function(p, d, q) {
    if (self$stage %in% c("estimate", "refit")) {
      p <- p[p <= floor(NROW(self$data) / 3)]
      q <- q[q <= floor(NROW(self$data) / 3)]
    }
    as.list(environment())
  },
  PDQ = function(P, D, Q, period = NULL) {
    stop("Seasonal VARIMA models are not yet supported.")
  },
  common_xregs,
  xreg = special_xreg(default_intercept = FALSE),
  .required_specials = c("pdq", "xreg"),
  .xreg_specials = names(common_xregs)
)

#' Estimate a VARIMA model
#'
#' Estimates a VARIMA model of a given order.
#'
#' Exogenous regressors and [`common_xregs`] can be specified in the model
#' formula.
#'
#' @aliases report.VARIMA
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param ... Further arguments for arima
#'
#' @section Specials:
#'
#' \subsection{pdq}{
#' The `pdq` special is used to specify non-seasonal components of the model.
#' \preformatted{
#' pdq(p = 0:5, d = 0:2, q = 0:5)
#' }
#'
#' \tabular{ll}{
#'   `p`      \tab The order of the non-seasonal auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#'   `d`      \tab The order of integration for non-seasonal differencing. If multiple values are provided, one of the values will be selected via repeated KPSS tests. \cr
#'   `q`      \tab The order of the non-seasonal moving average (MA) terms. If multiple values are provided, the one which minimises `ic` will be chosen. \cr
#' }
#' }

#' \subsection{xreg}{
#' Exogenous regressors can be included in an VARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
#'
#' The inclusion of a constant in the model follows the similar rules to [`stats::lm()`], where including `1` will add a constant and `0` or `-1` will remove the constant. If left out, the inclusion of a constant will be determined by minimising `ic`.
#'
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
#' [`MTS::VARMA()`], [`MTS:sVARMA()`].
#'
#' @examples
#'
#' eu_stocks <- as_tsibble(ts(EuStockMarkets, frequency = 1), pivot_longer = FALSE)
#'
#' eu_stocks |> 
#'   autoplot(vars(DAX, SMI, CAC, FTSE))
#'
#' fit <- eu_stocks %>%
#'   model(VARIMA(vars(DAX, SMI, CAC, FTSE) ~ pdq(1, 1, 1)))
#'
#' fit
#' @export
VARIMA <- function(formula, ...) {
  # ic <- match.arg(ic)
  # ic <- switch(ic, aicc = "AICc", aic = "AIC", bic = "BIC")
  varima_model <- new_model_class("VARIMA",
                                 train = train_varima,
                                 specials = specials_varima,
                                 origin = NULL,
                                 check = all_tsbl_checks
  )
  new_model_definition(varima_model, !!enquo(formula), ...)
}

#' @rdname VARIMA
#' @inheritParams forecast.VAR
#' @examples
#' 
#' fit %>%
#'   forecast(h = 50) %>%
#'   autoplot(tail(eu_stocks, 100))
#' @export
forecast.VARIMA <- function(object, new_data = NULL, specials = NULL,
                         bootstrap = FALSE, times = 5000, ...) {
  mts_varma_pred(object, NROW(new_data))
}

# nocov start
# Adapted from MTS::VARMApred to return distributional distributions
mts_varma_pred <- function(model, h) {
  x = as.matrix(model$data)
  resi = as.matrix(model$residuals)
  sig = model$Sigma
  Phi = model$Phi
  Theta = model$Theta
  Ph0 = model$Ph0
  p = model$ARorder
  q = model$MAorder
  if (p < 0) 
    p = 0
  if (q < 0) 
    q = 0
  if (h < 1) 
    h = 1
  nT = dim(x)[1]
  k = dim(x)[2]
  T1 = dim(resi)[1]
  if (nT > T1) {
    r1 = matrix(0, (nT - T1), k)
    resi = rbind(r1, resi)
  }
  if (length(Ph0) < 1) 
    Ph0 = rep(0, k)
  orig = nT
  px = x[1:orig, ]
  presi = resi[1:orig, ]
  psi = diag(rep(1, k))
  wk = c(psi)
  lag = max(1, h)
  for (i in 1:lag) {
    if (i <= p) {
      idx = (i - 1) * k
      tmp = Phi[, (idx + 1):(idx + k)]
    }
    else {
      tmp = matrix(0, k, k)
    }
    if (i <= q) {
      mdx = (i - 1) * k
      tmp = tmp - Theta[, (mdx + 1):(mdx + k)]
    }
    jj = i - 1
    jp = min(jj, p)
    if (jp > 0) {
      for (j in 1:jp) {
        jdx = (j - 1) * k
        idx = (i - j) * k
        w1 = Phi[, (jdx + 1):(jdx + k)]
        w2 = psi[, (idx + 1):(idx + k)]
        tmp = tmp + w1 %*% w2
      }
    }
    psi = cbind(psi, tmp)
    wk = cbind(wk, c(tmp))
  }
  mu = sigma = vector("list", h)
  for (j in 1:h) {
    fcst = Ph0
    Sig = sig
    t = orig + j
    if (p > 0) {
      for (ii in 1:p) {
        idx = (ii - 1) * k
        ph = Phi[, (idx + 1):(idx + k)]
        fcst = fcst + matrix(px[(t - ii), ], 1, k) %*% 
          t(ph)
      }
    }
    if (q > 0) {
      for (jj in 1:q) {
        idx = (jj - 1) * k
        if ((t - jj) <= orig) {
          th = Theta[, (idx + 1):(idx + k)]
          fcst = fcst - matrix(resi[(t - jj), ], 1, 
                               k) %*% t(th)
        }
      }
    }
    px = rbind(px, fcst)
    if (j > 1) {
      Sig = sig
      for (jj in 2:j) {
        jdx = (jj - 1) * k
        wk = psi[, (jdx + 1):(jdx + k)]
        Sig = Sig + wk %*% sig %*% t(wk)
      }
    }
    
    mu[[j]] <- fcst
    sigma[[j]] <- Sig
  }
  
  # Adjust to be forecasts in levels
  d <- NROW(model$y_end)
  if (d > 0) {
    mu <- do.call(rbind, mu)
    mu <- diffinv(mu, differences = d, xi = model$y_end)[-1,,drop=FALSE]
    mu <- unname(split(mu, row(mu)))
  }
  
  distributional::dist_multivariate_normal(mu, sigma)
}
# nocov end


#' @rdname VARIMA
#' @inheritParams fitted.VAR
#'
#' @examples
#'
#' fitted(fit)
#' @export
fitted.VARIMA <- function(object, ...) {
  d <- NROW(object$y_end)
  fit <- object$data - residuals(object)[seq_len(nrow(object$data)) + d,]
  
  if (d > 0) {
    fit[seq_len(object$ARorder),] <- object$data[seq_len(object$ARorder),]
    fit <- diffinv(fit, differences = d, xi = object$y_start)
    fit[seq_len(object$ARorder + d),] <- NA
  }
  fit
}

#' @rdname VARIMA
#' @inheritParams residuals.VAR
#' 
#' @examples
#' residuals(fit)
#' @export
residuals.VARIMA <- function(object, ...) {
  rbind(
    matrix(NA, nrow(object$data) - nrow(object$residuals) + nrow(object$y_end), ncol(object$residuals)),
    object$residuals
  )
}

#' @export
model_sum.VARIMA <- function(x) {
  sprintf(
    "VARIMA(%i,%i,%i)%s", 
    x$ARorder, NROW(x$y_end), x$MAorder, ifelse(x$cnst, " w/ mean", ""))
}

#' @rdname VARIMA
#' @inheritParams tidy.VAR
#'
#' @examples
#' tidy(fit)
#' @export
tidy.VARIMA <- function(x, ...) {
  colnames(x$coef) <- cn <- colnames(x$data)
  k <- NCOL(x$coef)
  rownames(x$coef) <- c(
    if(x$cnst) "constant" else NULL,
    paste0("lag(", rep(cn, x$ARorder), ",", rep(seq_len(x$ARorder), each = k), ")"),
    paste0("lag(e_", rep(cn, x$MAorder), ",", rep(seq_len(x$MAorder), each = k), ")")
  )
  
  rdf <- (x$ARorder + x$MAorder) * k^2 + k*(k+1)/2
  
  coef <- dplyr::as_tibble(x$coef, rownames = "term")
  coef <- tidyr::gather(coef, ".response", "estimate", !!!syms(colnames(x$coef)))
  mutate(
    coef,
    std.error = as.numeric(x$secoef),
    statistic = !!sym("estimate") / !!sym("std.error"),
    p.value = 2 * stats::pt(abs(!!sym("statistic")), rdf, lower.tail = FALSE)
  )
}


#' @rdname VARIMA
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @examples
#' glance(fit)
#' @export
glance.VARIMA <- function(x, ...) {
  tibble(
    sigma2 = list(x$Sigma),
    AIC = x$aic, BIC = x$bic
  )
}

#' @rdname VARIMA
#' 
#' @examples
#' report(fit)
#' 
#' @export
report.VARIMA <- function(object, ...) {
  coef <- tidy(object)
  coef <- map(
    split(coef, factor(coef$.response, levels = unique(coef$.response))),
    function(x) `colnames<-`(rbind(x$estimate, s.e. = x$std.error), x$term)
  )
  
  imap(coef, function(par, nm) {
    cat(sprintf("\nCoefficients for %s:\n", nm))
    print.default(round(par, digits = 4), print.gap = 2)
  })
  
  cat("\nResidual covariance matrix:\n")
  print.default(round(object$Sigma, 4))
  
  
  cat(
    sprintf(
      "\nAIC = %s\tBIC = %s",
      format(round(object$aic, 2L)),
      format(round(object$bic, 2L))
    )
  )
}

#' @rdname VARIMA
#' 
#' @examples
#' generate(fit, h = 10)
#' 
#' @export
generate.VARIMA <- function(x, new_data, specials, ...){
  if (!".innov" %in% names(new_data)) {
    new_data[[".innov"]] <- generate(
      distributional::dist_multivariate_normal(
        list(matrix(0, ncol = NCOL(x$Sigma))), list(x$Sigma)
      ),
      nrow(new_data)
    )[[1L]]
  }
  kr <- key_data(new_data)$.rows
  h <- lengths(kr)
  
  new_data$.sim <- do.call(
    rbind, lapply(kr, function(i) {
      mts_varma_sim(x, new_data[[".innov"]][i,,drop=FALSE])
    })
  )
  new_data
}


# nocov start
# Adapted from MTS::VARMAsim to simulate from model
mts_varma_sim <- function(model, innov) {
  cnst = if(model$cnst) model$Ph0 else 0
  phi = model$Phi
  nar = model$ARorder
  theta = model$Theta
  nma = model$MAorder
  k = NCOL(model$data)
  lastobs <- model$data[NROW(model$data) - seq_len(nar) + 1,,drop=FALSE]
  lastres <- model$residuals[NROW(model$residuals) - seq_len(nma) + 1,,drop=FALSE]
  zt <- matrix(numeric(length(innov)), ncol = ncol(innov))
  for (it in seq_len(NROW(innov))) {
    tmp = innov[it,,drop=FALSE]
    if (nma > 0) {
      for (j in seq_len(nma)) {
        jdx = (j - 1) * k
        thej = theta[, (jdx + 1):(jdx + k)]
        tmp = tmp - lastres[j,,drop=FALSE] %*% t(thej)
      }
    }
    if (nar > 0) {
      for (i in 1:nar) {
        idx = (i - 1) * k
        phj = phi[, (idx + 1):(idx + k)]
        tmp = tmp + lastobs[i,,drop=FALSE] %*% t(phj)
      }
    }
    lastobs[seq_len(nrow(lastobs)-1) + 1,] <- lastobs[seq_len(nrow(lastobs)-1),]
    lastobs[1,] = zt[it, ] = cnst + tmp
    lastres[seq_len(nrow(lastres)-1) + 1,] <- lastres[seq_len(nrow(lastres)-1),]
    lastres[1,] <- innov[it,]
  }
  
  d <- NROW(model$y_end)
  if (d > 0) {
    zt <- diffinv(zt, differences = d, xi = model$y_end)[-1,,drop=FALSE]
  }
  
  zt
}
# nocov end

#' @rdname VARIMA
#'
#' @param x A fitted model.
#' @param impulse A character string specifying the name of the variable that is shocked (the impulse variable).
#' @param orthogonal If TRUE, orthogonalised impulse responses will be computed.
#'
#' @examples
#' IRF(fit, h = 10, impulse = "FTSE")
#' 
#' @export
IRF.VARIMA <- function(x, new_data, specials, impulse = NULL, orthogonal = FALSE, ...) {
  # Zero out end of data
  x$data[seq_along(x$data)] <- 0
  x$residuals[seq_along(x$residuals)] <- 0
  x$y_end[seq_along(x$y_end)] <- 0
  
  # Remove regressors
  x$cnst <- FALSE
  
  # Add shocks
  if (".impulse" %in% names(new_data)) { 
    names(new_data)[names(new_data) == ".impulse"] <- ".innov"
  } else {
    new_data$.innov <- matrix(0, nrow = nrow(new_data), ncol = ncol(x$data),
                              dimnames = dimnames(x$data))
    new_data$.innov[1, impulse] <- 1
  }
  
  # Orthogonalised shocks
  if(orthogonal) {
    # Use Cholesky decomposition to orthogonalise the shocks / innovations
    new_data$.innov <- new_data$.innov %*% chol(x$Sigma)
  }
  
  irf <- generate(x, new_data, specials)
  irf[colnames(x$data)] <- split(irf$.sim, col(irf$.sim))
  irf$.innov <- irf$.sim <- NULL
  irf
}