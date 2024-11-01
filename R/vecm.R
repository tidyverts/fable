#' @importFrom stats ts
train_vecm <- function(.data, specials, ic, ...) {
  # Get args
  p <- specials$AR[[1]]$p

  # Get response variables
  y <- invoke(cbind, unclass(.data)[measured_vars(.data)])
  
  # Get xreg
  constant <- specials$xreg[[1]]$constant
  xreg <- Reduce(cbind, specials$xreg[[1]]$xreg) # Convert from df to matrix
  
  # Choose best model
  reduce(transpose(expand.grid(p = p, constant = constant)),
         function(best, args) {
           new <- estimate_vecm(y, args$p, xreg, args$constant, ...)
           if ((new$fit[[ic]] %||% Inf) < (best$fit[[ic]] %||% Inf)) {
             best <- new
           }
           best
         },
         .init = NULL
  )
}

estimate_vecm <- function(y, p, sr_xreg, constant, r, ...) {
  if (constant) {
    sr_xreg <- cbind(constant = rep(1, NROW(y)), sr_xreg)
  }
  
  dy <- diff(y)
  # if (p > 0) {
    # y <- y[-seq_len(p), , drop = FALSE]
    # dy <- dy[-seq_len(p), , drop = FALSE]
  # }
  sr_xreg <- sr_xreg[-seq_len(p + 1), , drop = FALSE]
  
  y_embed <- stats::embed(y, dimension = p + 1)
  dy_embed <- stats::embed(dy, dimension = p + 1)
  dy_y <- dy_embed[, seq_len(NCOL(dy)), drop = FALSE]
  y_y <- y_embed[-NROW(y_embed), seq_len(NCOL(y)), drop = FALSE]
  dy_x_sr <- dy_embed[, -seq_len(NCOL(dy)), drop = FALSE]
  
  colnames(dy_x_sr) <- pmap_chr(
    expand.grid(colnames(y), seq_len(p)), 
    sprintf, fmt = "lag(%s,%i)"
  )
  dm_sr <- cbind(dy_x_sr, sr_xreg)
  j <- complete.cases(dm_sr, dy_y, y_y)
  
  # Estimate covariance matrices with regression
  dy_fit <- stats::lm.fit(dm_sr[j, , drop = FALSE], dy_y[j, , drop = FALSE])
  y_fit <- stats::lm.fit(dm_sr[j, , drop = FALSE], y_y[j, , drop = FALSE])
  u <- residuals(dy_fit)
  v <- residuals(y_fit)
  
  S00 <- crossprod(u)
  S11 <- crossprod(v)
  S01 <- crossprod(u, v)
  
  # Solve the Eigenvalue problem
  # S_1^{−1} S_0 β= λ β
  eig_mat <- solve(S11) %*% t(S01) %*% solve(S00) %*% S01
  eig_dcmp <- eigen(eig_mat)
  eig_vec <- Re(eig_dcmp$vectors)
  eig_val <- Re(eig_dcmp$values)
  
  # Normalise the eigenvector matrix using Phillips triangular representation
  eig_vec_norm <- apply(eig_vec, 2, function(x) x/c(sqrt(t(x) %*% S11 %*% x)))
  # Rescale by columns to make diagonal elements 1
  eig_vec_norm <- t(t(eig_vec_norm)/diag(eig_vec_norm))
  
  # Compute the error correction term (beta' * Y_lagged)
  beta <- eig_vec_norm[, seq_len(r), drop = FALSE]
  C <- rbind(diag(r), matrix(0, nrow = nrow(eig_vec_norm) - r, ncol = r))
  beta <- beta %*% solve(t(C) %*% beta)
  beta[seq_len(r), seq_len(r)] <- diag(r)
  dimnames(beta) <- list(colnames(y), paste0("r", seq_len(r)))
  
  ECT <- y_y %*% beta
  colnames(ECT) <- paste0("ECT", seq_len(r))
  
  # Regress ECT and short run terms on Y
  fit <- stats::lm.fit(cbind(ECT, dm_sr[j, , drop = FALSE]), dy_y[j, , drop = FALSE])
  
  # Estimate the Adjustment Coefficients (Alpha)
  alpha <- solve(t(ECT) %*% ECT) %*% t(ECT) %*% dm_sr
  
  # Structure results
  resid <- matrix(NA_real_, nrow = nrow(y), ncol = ncol(y))
  resid[c(rep.int(FALSE, p + 1), j), ] <- fit$residuals
  if (is_empty(fit$coefficients)) {
    coef <- matrix(nrow = 0, ncol = NCOL(y))
  }
  else {
    coef <- as.matrix(fit$coefficients)
  }
  colnames(coef) <- colnames(y)
  
  nr <- sum(j)
  nc <- NCOL(y)
  sig <- crossprod(fit$residuals)
  sig_det <- det(sig / nr)
  loglik <- -(nr * nc / 2) * log(2 * pi) - (nr / 2) * log(sig_det) -
    (1 / 2) * sum(diag(resid %*% solve(sig / nr) %*% t(resid)), na.rm = TRUE)
  npar <- length(fit$coef[-seq_len(r),]) + 2 * nc * r - r^2
  aic <- -2 * loglik + 2 * npar
  bic <- aic + npar * (log(nr) - 2)
  aicc <- aic + 2 * npar * (npar + 1) / (nr - npar - 1)
  
  # Output model
  structure(
    list(
      coef = coef,
      fits = rbind(matrix(nrow = p, ncol = NCOL(y)), fit$fitted.values),
      resid = rbind(matrix(nrow = p, ncol = NCOL(y)), resid),
      fit = tibble(
        sigma2 = list(sig / fit$df.residual),
        alpha = list(alpha), beta = list(beta),
        log_lik = loglik, AIC = aic, AICc = aicc, BIC = bic
      ),
      spec = tibble(p = p, r = r, constant = constant),
      last_obs = y[NROW(y) - seq_len(p + 1) + 1, , drop = FALSE],
      model = fit
    ),
    class = "VECM"
  )
}

specials_vecm <- new_specials(
  AR = function(p = 0:5) {
    if (any(p < 0)) {
      warn("The AR order must be non-negative. Only non-negative orders will be considered.")
      p <- p[p >= 0]
    }
    list(p = p)
  },
  common_xregs,
  xreg = function(...) {
    dots <- enexprs(...)
    env <- map(enquos(...), get_env)
    env[map_lgl(env, compose(is_empty, env_parents))] <- NULL
    env <- if (!is_empty(env)) get_env(env[[1]]) else base_env()
    
    # Mask user defined lag to retain history when forecasting
    env <- env_bury(env, lag = lag)
    
    constants <- map_lgl(dots, inherits, "numeric")
    constant_given <- any(map_lgl(dots[constants], `%in%`, -1:1))
    
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    xreg <- model.frame(model_formula, data = env, na.action = stats::na.pass)
    list(
      constant = if (constant_given) as.logical(terms(xreg) %@% "intercept") else c(TRUE, FALSE),
      xreg = if (NCOL(xreg) == 0) NULL else xreg
    )
  },
  .required_specials = c("AR", "xreg"),
  .xreg_specials = names(common_xregs)
)

#' Estimate a VECM model
#'
#' Searches through the vector of lag orders to find the best VECM model which
#' has lowest AIC, AICc or BIC value. The model is estimated using the Johansen
#' procedure (maximum likelihood).
#'
#' Exogenous regressors and [`common_xregs`] can be specified in the model
#' formula.
#'
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param r The number of cointegrating relationships
#' @param ... Further arguments for arima
#'
#' @section Specials:
#'
#' \subsection{AR}{
#' The `AR` special is used to specify the lag order for the auto-regression.
#' \preformatted{
#' AR(p = 0:5)
#' }
#'
#' \tabular{ll}{
#'   `p`        \tab The order of the auto-regressive (AR) terms. If multiple values are provided, the one which minimises `ic` will be chosen.\cr
#' }
#' }
#'
#' \subsection{xreg}{
#' Exogenous regressors can be included in an VECM model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
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
#' @examples
#'
#' lung_deaths <- cbind(mdeaths, fdeaths) %>%
#'   as_tsibble(pivot_longer = FALSE)
#'
#' fit <- lung_deaths %>%
#'   model(VECM(vars(mdeaths, fdeaths) ~ AR(3)))
#'
#' report(fit)
#'
#' fit %>%
#'   forecast() %>%
#'   autoplot(lung_deaths)
#' @export
VECM <- function(formula, ic = c("aicc", "aic", "bic"), r = 1L, ...) {
  ic <- match.arg(ic)
  ic <- switch(ic, aicc = "AICc", aic = "AIC", bic = "BIC")
  vecm_model <- new_model_class("VECM",
                                 train = train_vecm,
                                 specials = specials_vecm,
                                 origin = NULL,
                                 check = all_tsbl_checks
  )
  new_model_definition(vecm_model, !!enquo(formula), ic = ic, r = r, ...)
}

#' @export
forecast.VECM <- function(object, new_data = NULL, specials = NULL,
                         bootstrap = FALSE, times = 5000, ...) {
  # Convert to VAR
  coef <- object$coef
  r <- object$spec$r
  p <- object$spec$p
  Pi <- object$fit$beta[[1]] %*% coef[seq_len(r),,drop = FALSE]
  k <- NCOL(coef)
  NROW(coef)
  var_coef <- matrix(0, nrow = NROW(coef) - r + k, ncol = k, )
  
  # Adjust AR coefficients
  coef_ar <- rbind(coef[seq(r + 1, r + p * k),,drop = FALSE], matrix(0, k, k))
  for (l in seq_len(p)) {
    var_coef[l*k + seq_len(k),] <- -(coef_ar[(l-1)*k + seq_len(k),] - coef_ar[l*k + seq_len(k),])
    var_coef[seq_len(k),] <- var_coef[seq_len(k),] + var_coef[l*k + seq_len(k),]
  }
  var_coef[seq_len(k),] <- Pi + diag(k) - var_coef[seq_len(k),]
  
  # Add regressors
  if ((NROW(coef) - r) > p * k) {
    var_coef[seq((p + 1) * k + 1, NROW(var_coef)),] <- coef[seq(p * k + r + 1, NROW(coef)),,drop = FALSE]
  }
  
  # Update model with VAR terms
  object$spec$p <- p + 1
  object$coef <- var_coef
  
  forecast.VAR(object, new_data, specials, ...)
}

#' @export
fitted.VECM <- function(object, ...) {
  object$fits
}

#' @export
residuals.VECM <- function(object, ...) {
  object$resid
}

#' @export
model_sum.VECM <- function(x) {
  sprintf("VECM(%s, r=%i)%s", x$spec$p, x$spec$r, ifelse(x$spec$constant, " w/ mean", ""))
}

#' @export
tidy.VECM <- function(x, ...) {
  rdf <- x$model$df.residual
  res <- split(x$resid, col(x$resid))
  rss <- map_dbl(res, function(resid) sum(resid^2, na.rm = TRUE))
  resvar <- rss / rdf
  rank <- x$model$rank
  
  R <- chol2inv(x$model$qr$qr[seq_len(rank), seq_len(rank), drop = FALSE])
  se <- map(resvar, function(resvar) sqrt(diag(R) * resvar))
  
  coef <- dplyr::as_tibble(x$coef, rownames = "term")
  coef <- tidyr::gather(coef, ".response", "estimate", !!!syms(colnames(x$coef)))
  mutate(
    coef,
    std.error = unlist(se),
    statistic = !!sym("estimate") / !!sym("std.error"),
    p.value = 2 * stats::pt(abs(!!sym("statistic")), rdf, lower.tail = FALSE)
  )
}


#' Glance a VECM
#'
#' Construct a single row summary of the VECM model.
#'
#' Contains the variance of residuals (`sigma2`), the log-likelihood
#' (`log_lik`), the cointegrating vector (`beta`) and information criterion
#' (`AIC`, `AICc`, `BIC`).
#'
#' @inheritParams generics::glance
#'
#' @return A one row tibble summarising the model's fit.
#'
#' @export
glance.VECM <- function(x, ...) {
  x$fit
}

#' @export
report.VECM <- function(object, ...) {
  coef <- tidy(object)
  coef <- map(
    split(coef, factor(coef$.response, levels = unique(coef$.response))),
    function(x) `colnames<-`(rbind(x$estimate, s.e. = x$std.error), x$term)
  )
  
  cat("\n Cointegrating vector:\n")
  print.default(round(object$fit$beta[[1]], 4))
  
  imap(coef, function(par, nm) {
    cat(sprintf("\nCoefficients for %s:\n", nm))
    print.default(round(par, digits = 4), print.gap = 2)
  })
  
  cat("\nResidual covariance matrix:\n")
  print.default(round(object$fit$sigma2[[1]], 4))
  
  cat(sprintf("\nlog likelihood = %s\n", format(round(object$fit$log_lik, 2L))))
  
  cat(
    sprintf(
      "AIC = %s\tAICc = %s\tBIC = %s",
      format(round(object$fit$AIC, 2L)),
      format(round(object$fit$AICc, 2L)),
      format(round(object$fit$BIC, 2L))
    )
  )
}

#' @inherit generate.ETS
#' 
#' @export
generate.VECM <- function(x, new_data, specials, ...){
  coef <- x$coef
  K <- NCOL(coef)
  if (!".innov" %in% names(new_data)) {
    new_data[[".innov"]] <- generate(distributional::dist_multivariate_normal(list(matrix(0, ncol = K)), x$fit$sigma2), nrow(new_data))[[1L]]
  }
  p <- x$spec$p
  kr <- key_data(new_data)$.rows
  h <- lengths(kr)
  
  # VECM params
  r <- x$spec$r
  alpha <-x$fit$alpha[[1]]
  beta <- x$fit$beta[[1]]
  
  # Get xreg
  xreg <- specials$xreg[[1]]$xreg
  
  # Generate paths
  vecm_sim <- function(i) {
    if (x$spec$constant) {
      xreg <- cbind(constant = rep_len(1, length(i)), xreg)
    }
    
    .innov <- new_data$.innov[i,]
    
    .sim <- matrix(NA, nrow = h, ncol = K)
    y_lag <- x$last_obs
    for (i in seq_along(i)) {
      # Difference (without diff() since matrix is in reverse order)
      Z <- t(y_lag[-NROW(y_lag),] - y_lag[-1,])
      if (!is.null(xreg)) {
        Z <- c(Z, t(xreg[i, ]))
      }
      
      # Error correction term
      ect <- t(y_lag[1,]) %*% beta %*% coef[seq_len(r),]
      
      # Short-run dynamics
      st <- c(Z) %*% coef[-seq_len(r),]
      
      .sim[i, ] <- y_lag[1,] + ect + st + .innov[i,]
      y_lag <- rbind(.sim[i, , drop = FALSE], y_lag)[seq_len(p + 1), , drop = FALSE]
    }
    
    .sim
  }
  
  new_data$.sim <- do.call(rbind, lapply(kr, vecm_sim))
  
  new_data
}

#' Calculate impulse responses from a fable model
#'
#' Simulates future paths from a dataset using a fitted model. Innovations are
#' sampled by the model's assumed error distribution. If `bootstrap` is `TRUE`,
#' innovations will be sampled from the model's residuals. If `new_data`
#' contains the `.innov` column, those values will be treated as innovations.
#'
#' @inheritParams forecast.ETS
#' @param x A fitted model.
#' @param impulse A character string specifying the name of the variable that is shocked (the impulse variable).
#' @param orthogonal If TRUE, orthogonalised impulse responses will be computed.
#'
#' @export
IRF.VECM <- function(x, new_data, specials, impulse = NULL, orthogonal = FALSE, ...) {
  # Zero out end of data
  x$last_obs[seq_along(x$last_obs)] <- 0
  
  # Remove regressors
  n_vecm <- x$spec$r + x$spec$p * ncol(x$coef)
  if(nrow(x$coef) > n_vecm) {
    x$coef[seq(n_vecm + 1, nrow(x$coef)),] <- 0
  }
  
  # Add shocks
  if (".impulse" %in% names(new_data)) { 
    names(new_data)[names(new_data) == ".impulse"] <- ".innov"
  } else {
    new_data$.innov <- matrix(0, nrow = nrow(new_data), ncol = ncol(x$last_obs),
                              dimnames = dimnames(x$last_obs))
    new_data$.innov[1, impulse] <- 1
  }
  
  # Orthogonalised shocks
  if(orthogonal) {
    # Use Cholesky decomposition to orthogonalise the shocks / innovations
    new_data$.innov <- new_data$.innov %*% chol(x$fit$sigma2[[1L]])
  }
  
  irf <- generate(x, new_data, specials)
  irf[colnames(x$coef)] <- split(irf$.sim, col(irf$.sim))
  irf$.innov <- irf$.sim <- NULL
  irf
}