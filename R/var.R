#' @importFrom stats ts
train_var <- function(.data, formula, specials, ...){
  # Get args
  p <- specials$AR[[1]]$p
  
  # Get response variables
  y <- as.matrix(.data[measured_vars(.data)])
  
  # Get xreg
  xreg <- specials$xreg[[1]]
  
  # Compute response lags
  y_lag <- stats::embed(y, dimension = p + 1)[, -(seq_len(NCOL(y)))]
  colnames(y_lag) <- pmap_chr(expand.grid(colnames(y), seq_len(p)), 
                              sprintf, fmt = "lag(%s,%i)")
  if(p > 0){
    xreg <- xreg[-seq_len(p),, drop = FALSE]
    y <- y[-seq_len(p),, drop = FALSE]
  }
  dm <- cbind(y_lag, xreg)
  fit <- stats::lm.fit(as.matrix(dm), y)
  
  resid <- as.matrix(fit$residuals)
  
  nr <- NROW(y)
  nc <- NCOL(y)
  sig <- crossprod(fit$residuals)
  sig_det <- det(sig/nr)
  loglik <- -(nr * nc/2) * log(2 * pi) - (nr/2) * log(sig_det) - 
    (1/2) * sum(diag(resid %*% solve(sig/nr) %*% t(resid)))
  
  aic <- log(sig_det) + (2/nr) * (p * nc^2 + nc * NCOL(xreg))
  hq <- log(sig_det) + (2 * log(log(nr))/nr) * (p * nc^2 + nc * NCOL(xreg))
  sc <- log(sig_det) + (log(nr)/nr) * (p * nc^2 + nc * NCOL(xreg))
  fpe <- ((nr + NCOL(dm))/(nr - NCOL(dm)))^nc * sig_det
  
  # Output model
  structure(
    list(
      coef = as.matrix(fit$coefficients),
      fits = rbind(matrix(nrow = p, ncol = NCOL(y)), y - resid),
      resid = rbind(matrix(nrow = p, ncol = NCOL(y)), resid),
      fit = tibble(sigma2 = list(sig/fit$df.residual), logLik = loglik,
                   aic = aic, hq = hq, sc = sc, fpe = fpe),
      spec = tibble(p = p),
      last_obs = y[NROW(y) - seq_len(p) + 1,,drop = FALSE],
      model = fit
    ),
    class = "VAR"
  )
}

specials_var <- new_specials(
  AR = function(p = 0){
    list(p=p)
  },
  common_xregs,
  xreg = function(...){
    model_formula <- new_formula(
      lhs = NULL,
      rhs = reduce(enexprs(...), function(.x, .y) call2("+", .x, .y))
    )
    mf <- model.frame(model_formula, data = self$data, na.action = stats::na.pass)
    if(mf%@%"terms"%@%"intercept" == 1){
      mf <- cbind(constant = rep(1, NROW(self$data)), mf)
    }
    mf
  },
  .required_specials = c("AR", "xreg"),
  .xreg_specials = names(common_xregs)
)

#' Estimate an VAR model
#' 
#' @param formula Model specification (see "Specials" section).
#' @param ... Further arguments for arima
#' 
#' @section Specials:
#' 
#' \subsection{pdq}{
#' The `AR` special is used to specify the lag order for the auto-regression.
#' \preformatted{
#' AR(p = 0)
#' }
#'
#' \tabular{ll}{
#'   `p`        \tab The order of the auto-regressive (AR) terms.\cr
#' }
#' }
#' 
#' #' \subsection{xreg}{
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
#' @examples 
#' 
#' lung_deaths <- cbind(mdeaths, fdeaths) %>%
#'   as_tsibble(pivot_longer = FALSE)
#'   
#' fit <- lung_deaths %>%
#'   model(VAR(vars(log(mdeaths), fdeaths) ~ AR(3)))
#' fit
#'
#' fit %>% 
#'   forecast()
#' 
#' @export
VAR <- function(formula, ...){
  varma_model <- new_model_class("VAR", train = train_var, 
                                 specials = specials_var,
                                 origin = NULL,
                                 check = all_tsbl_checks)
  new_model_definition(varma_model, !!enquo(formula), ...)
}

#' @export
forecast.VAR <- function(object, new_data = NULL, specials = NULL, 
                           bootstrap = FALSE, times = 5000, ...){
  if(!is_regular(new_data)){
    abort("Forecasts must be regularly spaced.")
  }
  if(bootstrap){
    abort("Bootstrapped forecasts for VARs are not yet implemented.")
  }
  
  # Get xreg
  xreg <- specials$xreg[[1]]
  
  h <- NROW(new_data)
  p <- object$spec$p
  coef <- object$coef
  K <- NCOL(coef)
  
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
    Z <-  c(t(y_lag), xreg[i,])
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
  sprintf("VAR(%s)", x$spec$p)
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