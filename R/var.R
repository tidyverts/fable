#' @importFrom stats ts
train_var <- function(.data, formula, specials, ic, ...){
  # Get args
  p <- specials$AR[[1]]$p
  
  # Get response variables
  y <- as.matrix(.data[measured_vars(.data)])
  
  # Get xreg
  constant <- specials$xreg[[1]]$constant
  xreg <- specials$xreg[[1]]$xreg
  
  # Choose best model
  reduce(transpose(expand.grid(p = p, constant = constant)),
         function(best, args){
           new <- estimate_var(y, args$p, xreg, args$constant)
           if((new$fit[[ic]] %||% Inf) < (best$fit[[ic]] %||% Inf)){
             best <- new
           }
           best
         }, .init = NULL)
}

estimate_var <- function(y, p, xreg, constant){
  if(constant){
    xreg <- cbind(constant = rep(1, NROW(y)), xreg)
  }
  
  y_lag <- stats::embed(y, dimension = p + 1)[, -(seq_len(NCOL(y))), drop = FALSE]
  colnames(y_lag) <- pmap_chr(expand.grid(colnames(y), seq_len(p)), 
                              sprintf, fmt = "lag(%s,%i)")
  if(p > 0){
    xreg <- xreg[-seq_len(p),, drop = FALSE]
    y <- y[-seq_len(p),, drop = FALSE]
  }
  dm <- cbind(y_lag, xreg)
  fit <- stats::lm.fit(as.matrix(dm), y)
  
  resid <- as.matrix(fit$residuals)
  
  nr <- NROW(stats::na.omit(y))
  nc <- NCOL(y)
  sig <- crossprod(fit$residuals)
  sig_det <- det(sig/nr)
  loglik <- -(nr * nc/2) * log(2 * pi) - (nr/2) * log(sig_det) - 
    (1/2) * sum(diag(resid %*% solve(sig/nr) %*% t(resid)))
  
  npar <- (length(fit$coef) + nc^2)
  aic <- -2*loglik + 2*npar
  bic <- aic + npar * (log(nr) - 2)
  aicc <- aic + 2 * npar * (npar + 1) / (nr - npar - 1)
  
  # Output model
  structure(
    list(
      coef = as.matrix(fit$coefficients),
      fits = rbind(matrix(nrow = p, ncol = NCOL(y)), y - resid),
      resid = rbind(matrix(nrow = p, ncol = NCOL(y)), resid),
      fit = tibble(sigma2 = list(sig/fit$df.residual), logLik = loglik,
                   aic = aic, aicc = aicc, bic = bic),
      spec = tibble(p = p, constant = constant),
      last_obs = y[NROW(y) - seq_len(p) + 1,,drop = FALSE],
      model = fit
    ),
    class = "VAR"
  )
}

specials_var <- new_specials(
  AR = function(p = 0:5){
    if(any(p < 0)){
      warn("The AR order must be non-negative. Only non-negative orders will be considered.")
      p <- p[p >= 0]
    }
    list(p=p)
  },
  common_xregs,
  xreg = function(...){
    dots <- enexprs(...)
    
    constants <- map_lgl(dots, inherits, "numeric") 
    constant_given <- any(map_lgl(dots[constants], `%in%`, -1:1))
    
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(dots, function(.x, .y) call2("+", .x, .y))
    )
    xreg <- model.frame(model_formula, data = self$data, na.action = stats::na.pass)
    
    list(
      constant = if(constant_given) as_logical(terms(xreg)%@%"intercept") else c(TRUE, FALSE),
      xreg = if(NCOL(xreg) == 0) NULL else xreg
    )
  },
  .required_specials = c("AR", "xreg"),
  .xreg_specials = names(common_xregs)
)

#' Estimate an VAR model
#' 
#' @param formula Model specification (see "Specials" section).
#' @param ic The information criterion used in selecting the model.
#' @param ... Further arguments for arima
#' 
#' @section Specials:
#' 
#' \subsection{pdq}{
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
#' Exogenous regressors can be included in an ARIMA model without explicitly using the `xreg()` special. Common exogenous regressor specials as specified in [`common_xregs`] can also be used. These regressors are handled using [stats::model.frame()], and so interactions and other functionality behaves similarly to [stats::lm()].
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
#' @examples 
#' 
#' lung_deaths <- cbind(mdeaths, fdeaths) %>%
#'   as_tsibble(pivot_longer = FALSE)
#'   
#' fit <- lung_deaths %>%
#'   model(VAR(vars(log(mdeaths), fdeaths) ~ AR(3)))
#' 
#' report(fit)
#'
#' fit %>% 
#'   forecast() %>% 
#'   autoplot(lung_deaths)
#' 
#' @export
VAR <- function(formula, ic = c("aicc", "aic", "bic"), ...){
  ic <- match.arg(ic)
  varma_model <- new_model_class("VAR", train = train_var, 
                                 specials = specials_var,
                                 origin = NULL,
                                 check = all_tsbl_checks)
  new_model_definition(varma_model, !!enquo(formula), ic = ic, ...)
}

#' @export
forecast.VAR <- function(object, new_data = NULL, specials = NULL, 
                           bootstrap = FALSE, times = 5000, ...){
  if(bootstrap){
    abort("Bootstrapped forecasts for VARs are not yet implemented.")
  }
  
  h <- NROW(new_data)
  p <- object$spec$p
  coef <- object$coef
  K <- NCOL(coef)
  # Get xreg
  xreg <- specials$xreg[[1]]$xreg
  if(object$spec$constant){
    xreg <- cbind(constant = rep(1, h), xreg)
  }
  
  # Compute phi
  As <- rep(list(matrix(0, nrow = K, ncol = K)), max(h, p))
  for(i in seq_len(p)){
    As[[i]] <- coef[seq_len(K) + (i-1)*K,]
  }
  phi <- rep(list(matrix(0, nrow = K, ncol = K)), h + 1)
  phi[[1]] <- diag(K)
  for(i in seq_len(h) + 1){
    tmp1 <- phi[[1]] %*% As[[i - 1]]
    tmp2 <- matrix(0, nrow = K, ncol = K)
    idx <- rev(seq_len(i-2))
    for (j in seq_len(i-2)) {
      tmp2 <- tmp2 + phi[[j + 1]] %*% As[[idx[j]]]
    }
    phi[[i]] <- tmp1 + tmp2
  }
  # Compute sigma
  sigma.u <- object$fit$sigma2[[1]]
  sigma <- rep(list(matrix(nrow = K, ncol = K)), h)
  sigma[[1]] <- sigma.u
  
  for (i in seq_len(h-1) + 1) {
    adjust <- matrix(0, nrow = K, ncol = K)
    for (j in 2:i) {
      adjust <- adjust + t(phi[[j]]) %*% sigma.u %*% phi[[j]]
    }
    sigma[[i]] <- adjust + sigma[[1]]
  }
  sigma <- map(sigma, sqrt)
  
  # Compute forecasts
  fc <- matrix(NA, ncol = K, nrow = h)
  
  y_lag <- object$last_obs
  for(i in seq_len(h)){
    Z <-  c(t(y_lag), t(xreg[i,]))
    fc[i,] <- t(coef) %*% Z
    y_lag <- rbind(fc[i,,drop=FALSE], y_lag)[seq_len(p),,drop=FALSE]
  }
  
  if(NCOL(fc) == 1){
    dist <- dist_normal(fc, map_dbl(sigma, `[`, 1, 1))
  }
  else {
    dist <- dist_mv_normal(split(fc, row(fc)), sigma)
  }
  
  # Output forecasts
  construct_fc(
    split(fc, col(fc)), 
    map(seq_len(K), function(x,y) map_dbl(y, `[[`, x), map(sigma, diag)),
    dist
  )
}

#' @export
fitted.VAR <- function(object, ...){
  object$fits
}

#' @export
residuals.VAR <- function(object, ...){
  object$resid
}

#' @export
model_sum.VAR <- function(x){
  sprintf("VAR(%s)%s", x$spec$p, ifelse(x$spec$constant, " w/ mean", ""))
}

#' @export
tidy.VAR <- function(x){
  rdf <- x$model$df.residual
  res <- split(x$resid, col(x$resid))
  rss <- map_dbl(res, function(resid) sum(resid^2, na.rm = TRUE))
  resvar <- rss/rdf
  rank <- x$model$rank
  
  R <- chol2inv(x$model$qr$qr[seq_len(rank), seq_len(rank), drop = FALSE])
  se <- map(resvar, function(resvar) sqrt(diag(R) * resvar))
  
  dplyr::as_tibble(x$coef, rownames = "term") %>% 
    tidyr::gather(".response", "estimate", !!!syms(colnames(x$coef))) %>% 
    mutate(std.error = unlist(se),
           statistic = !!sym("estimate") / !!sym("std.error"),
           p.value = 2 * stats::pt(abs(!!sym("statistic")), rdf, lower.tail = FALSE))
}

#' @export
glance.VAR <- function(x, ...){
  x$fit
}

#' @export
report.VAR <- function(object, ...){
  coef <- tidy(object)
  coef <- map(
    split(coef, factor(coef$.response, levels = unique(coef$.response))),
    function(x) `colnames<-`(rbind(x$estimate, s.e. = x$std.error), x$term)
  )
  
  imap(coef, function(par, nm){
    cat(sprintf("\nCoefficients for %s:\n", nm))
    print.default(round(par, digits = 4), print.gap = 2)
  })
  
  cat("\nResidual covariance matrix:\n")
  print.default(round(object$fit$sigma2[[1]], 4))
  
  cat(sprintf("\nlog likelihood = %s\n", format(round(object$fit$logLik, 2L))))
  
  cat(
    sprintf("AIC = %s\tAICc = %s\tBIC = %s",
            format(round(object$fit$aic, 2L)),
            format(round(object$fit$aicc, 2L)),
            format(round(object$fit$bic, 2L))
    )
  )
}