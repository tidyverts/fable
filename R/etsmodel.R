etsmodel <- function(y, m, errortype, trendtype, seasontype, damped,
                     alpha=NULL, beta=NULL, gamma=NULL, phi=NULL,
                     lower, upper, opt.crit, nmse, bounds, maxit=2000, trace=FALSE) {
  if (seasontype ==  "N") {
    m <- 1
  }
  
  # Modify limits if alpha, beta or gamma have been specified.
  if (!is.null(alpha)) {
    upper[2] <- min(alpha, upper[2])
    upper[3] <- min(1 - alpha, upper[3])
  }
  if (!is.null(beta)) {
    lower[1] <- max(beta, lower[1])
  }
  if (!is.null(gamma)) {
    upper[1] <- min(1 - gamma, upper[1])
  }
  
  # Initialize smoothing parameters
  par <- initparam(alpha, beta, gamma, phi, trendtype, seasontype, damped, lower, upper, m)
  names(alpha) <- names(beta) <- names(gamma) <- names(phi) <- NULL
  par.noopt <- c(alpha = alpha, beta = beta, gamma = gamma, phi = phi)
  if (!is.null(par.noopt)) {
    par.noopt <- c(stats::na.omit(par.noopt))
  }
  if (!is.na(par["alpha"])) {
    alpha <- par["alpha"]
  }
  if (!is.na(par["beta"])) {
    beta <- par["beta"]
  }
  if (!is.na(par["gamma"])) {
    gamma <- par["gamma"]
  }
  if (!is.na(par["phi"])) {
    phi <- par["phi"]
  }
  
  if (!check.param(alpha, beta, gamma, phi, lower, upper, bounds, m)) {
    print(paste("Model: ETS(", errortype, ",", trendtype, ifelse(damped, "d", ""), ",", seasontype, ")", sep = ""))
    stop("Parameters out of range")
  }
  
  # Initialize state
  init.state <- initstate(y, m, trendtype, seasontype)
  nstate <- length(init.state)
  par <- c(par, init.state)
  lower <- c(lower, rep(-Inf, nstate))
  upper <- c(upper, rep(Inf, nstate))
  
  np <- length(par)
  if (np >= length(y) - 1) { # Not enough data to continue
    return(list(aic = Inf, bic = Inf, aicc = Inf, mse = Inf, amse = Inf, fit = NULL, par = par, states = init.state))
  }
  
  env <- etsTargetFunctionInit(
    par = par, y = y, nstate = nstate, errortype = errortype, trendtype = trendtype,
    seasontype = seasontype, damped = damped, par.noopt = par.noopt, lowerb = lower, upperb = upper,
    opt.crit = opt.crit, nmse = as.integer(nmse), bounds = bounds, m = m, pnames = names(par), pnames2 = names(par.noopt)
  )
  
  fred <- .Call(
    "etsNelderMead", par, env, -Inf,
    sqrt(.Machine$double.eps), 1.0, 0.5, 2.0, trace, maxit, PACKAGE = "fable"
  )
  
  fit.par <- fred$par
  
  names(fit.par) <- names(par)
  
  init.state <- fit.par[(np - nstate + 1):np]
  # Add extra state
  if (seasontype != "N") {
    init.state <- c(init.state, m * (seasontype == "M") - sum(init.state[(2 + (trendtype != "N")):nstate]))
  }
  
  if (!is.na(fit.par["alpha"])) {
    alpha <- fit.par["alpha"]
  }
  if (!is.na(fit.par["beta"])) {
    beta <- fit.par["beta"]
  }
  if (!is.na(fit.par["gamma"])) {
    gamma <- fit.par["gamma"]
  }
  if (!is.na(fit.par["phi"])) {
    phi <- fit.par["phi"]
  }
  
  e <- pegelsresid.C(y, m, init.state, errortype, trendtype, seasontype, damped, alpha, beta, gamma, phi, nmse)
  
  np <- np + 1
  ny <- length(y)
  aic <- e$lik + 2 * np
  bic <- e$lik + log(ny) * np
  aicc <- aic + 2 * np * (np + 1) / (ny - np - 1)
  
  mse <- e$amse[1]
  amse <- mean(e$amse)
  
  states <- e$states
  
  states_cn <- "l"
  if (trendtype != "N") {
    states_cn <- c(states_cn, "b")
  }
  if (seasontype != "N") {
    states_cn <- c(states_cn, paste("s", 1:m, sep = ""))
  }
  colnames(states) <- states_cn
  
  fit.par <- c(fit.par, par.noopt)
  if (errortype == "A") {
    fits <- y - e$e
  } else {
    fits <- y / (1 + e$e)
  }
  
  return(list(
    loglik = -0.5 * e$lik, aic = aic, bic = bic, aicc = aicc, mse = mse, amse = amse,
    fit = fred, residuals = e$e, fitted = fits,
    states = states, par = fit.par
  ))
}


etsTargetFunctionInit <- function(par, y, nstate, errortype, trendtype, seasontype, damped, par.noopt, lowerb, upperb,
                                  opt.crit, nmse, bounds, m, pnames, pnames2) {
  names(par) <- pnames
  names(par.noopt) <- pnames2
  alpha <- c(par["alpha"], par.noopt["alpha"])["alpha"]
  if (is.na(alpha)) {
    stop("alpha problem!")
  }
  if (trendtype != "N") {
    beta <- c(par["beta"], par.noopt["beta"])["beta"]
    if (is.na(beta)) {
      stop("beta Problem!")
    }
  }
  else {
    beta <- NULL
  }
  if (seasontype != "N") {
    gamma <- c(par["gamma"], par.noopt["gamma"])["gamma"]
    if (is.na(gamma)) {
      stop("gamma Problem!")
    }
  }
  else {
    m <- 1
    gamma <- NULL
  }
  if (damped) {
    phi <- c(par["phi"], par.noopt["phi"])["phi"]
    if (is.na(phi)) {
      stop("phi Problem!")
    }
  }
  else {
    phi <- NULL
  }
  
  # determine which values to optimize and which ones are given by the user/not needed
  optAlpha <- !is.null(alpha)
  optBeta <- !is.null(beta)
  optGamma <- !is.null(gamma)
  optPhi <- !is.null(phi)
  
  givenAlpha <- FALSE
  givenBeta <- FALSE
  givenGamma <- FALSE
  givenPhi <- FALSE
  
  if (!is.null(par.noopt["alpha"])) {
    if (!is.na(par.noopt["alpha"])) {
      optAlpha <- FALSE
      givenAlpha <- TRUE
    }
  }
  if (!is.null(par.noopt["beta"])) {
    if (!is.na(par.noopt["beta"])) {
      optBeta <- FALSE
      givenBeta <- TRUE
    }
  }
  if (!is.null(par.noopt["gamma"])) {
    if (!is.na(par.noopt["gamma"])) {
      optGamma <- FALSE
      givenGamma <- TRUE
    }
  }
  if (!is.null(par.noopt["phi"])) {
    if (!is.na(par.noopt["phi"])) {
      optPhi <- FALSE
      givenPhi <- TRUE
    }
  }
  
  
  if (!damped) {
    phi <- 1
  }
  if (trendtype == "N") {
    beta <- 0
  }
  if (seasontype == "N") {
    gamma <- 0
  }
  
  env <- new.env()
  
  res <- .Call(
    "etsTargetFunctionInit", y = y, nstate = nstate, errortype = switch(errortype, "A" = 1, "M" = 2),
    trendtype = switch(trendtype, "N" = 0, "A" = 1, "M" = 2), seasontype = switch(seasontype, "N" = 0, "A" = 1, "M" = 2),
    damped = damped, lowerb = lowerb, upperb = upperb,
    opt.crit = opt.crit, nmse = as.integer(nmse), bounds = bounds, m = m,
    optAlpha, optBeta, optGamma, optPhi,
    givenAlpha, givenBeta, givenGamma, givenPhi,
    alpha, beta, gamma, phi, env, PACKAGE = "fable"
  )
  res
}

initparam <- function(alpha, beta, gamma, phi, trendtype, seasontype, damped, lower, upper, m) {
  if (any(lower > upper)) {
    stop("Inconsistent parameter boundaries")
  }
  
  # Select alpha
  if (is.null(alpha)) {
    alpha <- lower[1] + 0.2 * (upper[1] - lower[1]) / m
    if (alpha > 1 || alpha < 0) {
      alpha <- lower[1] + 2e-3
    }
    par <- c(alpha = alpha)
  }
  else {
    par <- numeric(0)
  }
  
  # Select beta
  if (trendtype != "N" && is.null(beta)) {
    # Ensure beta < alpha
    upper[2] <- min(upper[2], alpha)
    beta <- lower[2] + 0.1 * (upper[2] - lower[2])
    if (beta < 0 || beta > alpha) {
      beta <- alpha - 1e-3
    }
    par <- c(par, beta = beta)
  }
  
  # Select gamma
  if (seasontype != "N" && is.null(gamma)) {
    # Ensure gamma < 1-alpha
    upper[3] <- min(upper[3], 1 - alpha)
    gamma <- lower[3] + 0.05 * (upper[3] - lower[3])
    if (gamma < 0 || gamma > 1 - alpha) {
      gamma <- 1 - alpha - 1e-3
    }
    par <- c(par, gamma = gamma)
  }
  
  # Select phi
  if (damped && is.null(phi)) {
    phi <- lower[4] + .99 * (upper[4] - lower[4])
    if (phi < 0 || phi > 1) {
      phi <- upper[4] - 1e-3
    }
    par <- c(par, phi = phi)
  }
  
  return(par)
}

check.param <- function(alpha, beta, gamma, phi, lower, upper, bounds, m) {
  if (bounds != "admissible") {
    if (!is.null(alpha)) {
      if (alpha < lower[1] || alpha > upper[1]) {
        return(0)
      }
    }
    if (!is.null(beta)) {
      if (beta < lower[2] || beta > alpha || beta > upper[2]) {
        return(0)
      }
    }
    if (!is.null(phi)) {
      if (phi < lower[4] || phi > upper[4]) {
        return(0)
      }
    }
    if (!is.null(gamma)) {
      if (gamma < lower[3] || gamma > 1 - alpha || gamma > upper[3]) {
        return(0)
      }
    }
  }
  if (bounds != "usual") {
    if (!admissible(alpha, beta, gamma, phi, m)) {
      return(0)
    }
  }
  return(1)
}

initstate <- function(y, m, trendtype, seasontype) {
  if (seasontype != "N") {
    # Do decomposition
    n <- length(y)
    if (n < 4) {
      stop("You've got to be joking (not enough data).")
    } else if (n < 3 * m) # Fit simple Fourier model.
    {
      fouriery <- fourier(seq_along(y), m, 1)
      trendy <- seq_along(y)
      fit <- stats::lm(y ~ trendy + fouriery)
      if (seasontype == "A") {
        y.d <- list(seasonal = y - fit$coef[1] - fit$coef[2] * (1:n))
      } else { # seasontype=="M". Biased method, but we only need a starting point
        y.d <- list(seasonal = y / (fit$coef[1] + fit$coef[2] * (1:n)))
      }
    }
    else { # n is large enough to do a decomposition
      y.d <- stats::decompose(stats::ts(y, frequency = m), type = switch(seasontype, A = "additive", M = "multiplicative"))
    }
    
    init.seas <- rev(y.d$seasonal[2:m]) # initial seasonal component
    names(init.seas) <- paste("s", 0:(m - 2), sep = "")
    # Seasonally adjusted data
    if (seasontype == "A") {
      y.sa <- y - as.numeric(y.d$seasonal)
    } else {
      init.seas <- pmax(init.seas, 1e-2) # We do not want negative seasonal indexes
      if (sum(init.seas) > m) {
        init.seas <- init.seas / sum(init.seas + 1e-2)
      }
      y.sa <- y / pmax(y.d$seasonal, 1e-2)
    }
  }
  else # non-seasonal model
  {
    m <- 1
    init.seas <- NULL
    y.sa <- y
  }
  
  maxn <- min(max(10, 2 * m), length(y.sa))
  if (trendtype == "N") {
    l0 <- mean(y.sa[1:maxn])
    b0 <- NULL
  }
  else # Simple linear regression on seasonally adjusted data
  {
    fit <- stats::lsfit(1:maxn, y.sa[1:maxn])
    if (trendtype == "A") {
      l0 <- fit$coef[1]
      b0 <- fit$coef[2]
      # If error type is "M", then we don't want l0+b0=0.
      # So perturb just in case.
      if (abs(l0 + b0) < 1e-8) {
        l0 <- l0 * (1 + 1e-3)
        b0 <- b0 * (1 - 1e-3)
      }
    }
    else # if(trendtype=="M")
    {
      l0 <- fit$coef[1] + fit$coef[2] # First fitted value
      if (abs(l0) < 1e-8) {
        l0 <- 1e-7
      }
      b0 <- (fit$coef[1] + 2 * fit$coef[2]) / l0 # Ratio of first two fitted values
      l0 <- l0 / b0 # First fitted value divided by b0
      if (abs(b0) > 1e10) { # Avoid infinite slopes
        b0 <- sign(b0) * 1e10
      }
      if (l0 < 1e-8 || b0 < 1e-8) # Simple linear approximation didn't work.
      {
        l0 <- max(y.sa[1], 1e-3)
        b0 <- max(y.sa[2] / y.sa[1], 1e-3)
      }
    }
  }
  
  names(l0) <- "l"
  if (!is.null(b0)) {
    names(b0) <- "b"
  }
  return(c(l0, b0, init.seas))
}

pegelsresid.C <- function(y, m, init.state, errortype, trendtype, seasontype, damped, alpha, beta, gamma, phi, nmse) {
  n <- length(y)
  p <- length(init.state)
  x <- numeric(p * (n + 1))
  x[1:p] <- init.state
  e <- numeric(n)
  lik <- 0
  if (!damped) {
    phi <- 1
  }
  if (trendtype == "N") {
    beta <- 0
  }
  if (seasontype == "N") {
    gamma <- 0
  }
  
  amse <- numeric(nmse)
  
  Cout <- .C(
    "etscalc",
    as.double(y),
    as.integer(n),
    as.double(x),
    as.integer(m),
    as.integer(switch(errortype, "A" = 1, "M" = 2)),
    as.integer(switch(trendtype, "N" = 0, "A" = 1, "M" = 2)),
    as.integer(switch(seasontype, "N" = 0, "A" = 1, "M" = 2)),
    as.double(alpha),
    as.double(beta),
    as.double(gamma),
    as.double(phi),
    as.double(e),
    as.double(lik),
    as.double(amse),
    as.integer(nmse),
    PACKAGE = "fable"
  )
  if (!is.na(Cout[[13]])) {
    if (abs(Cout[[13]] + 99999) < 1e-7) {
      Cout[[13]] <- NA
    }
  }
  
  return(list(lik = Cout[[13]], amse = Cout[[14]], e = Cout[[12]], states = matrix(Cout[[3]], nrow = n + 1, ncol = p, byrow = TRUE)))
}

admissible <- function(alpha, beta, gamma, phi, m) {
  if (is.null(phi)) {
    phi <- 1
  }
  if (phi < 0 || phi > 1 + 1e-8) {
    return(0)
  }
  if (is.null(gamma)) {
    if (alpha < 1 - 1 / phi || alpha > 1 + 1 / phi) {
      return(0)
    }
    if (!is.null(beta)) {
      if (beta < alpha * (phi - 1) || beta > (1 + phi) * (2 - alpha)) {
        return(0)
      }
    }
  }
  else if (m > 1) # Seasonal model
  {
    if (is.null(beta)) {
      beta <- 0
    }
    if (gamma < max(1 - 1 / phi - alpha, 0) || gamma > 1 + 1 / phi - alpha) {
      return(0)
    }
    if (alpha < 1 - 1 / phi - gamma * (1 - m + phi + phi * m) / (2 * phi * m)) {
      return(0)
    }
    if (beta < -(1 - phi) * (gamma / m + alpha)) {
      return(0)
    }
    
    # End of easy tests. Now use characteristic equation
    P <- c(phi * (1 - alpha - gamma), alpha + beta - alpha * phi + gamma - 1, rep(alpha + beta - alpha * phi, m - 2), (alpha + beta - phi), 1)
    roots <- polyroot(P)
    
    # cat("maxpolyroots: ", max(abs(roots)), "\n")
    
    if (max(abs(roots)) > 1 + 1e-10) {
      return(0)
    }
  }
  # Passed all tests
  return(1)
}