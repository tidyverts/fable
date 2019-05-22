#' @importFrom stats sd
train_mean <- function(.data, formula, specials, ...){
  y <- .data[[measured_vars(.data)]]
  
  if(all(is.na(y))){
    abort("All observations are missing, a model cannot be estimated without data.")
  }
  
  n <- length(y)
  y_mean <- mean(y, na.rm = TRUE)
  fits <- rep(y_mean, n)
  res <- y - fits
  sigma <- sd(y, na.rm = TRUE)
  
  structure(
    list(
      par = tibble(term = "mean", estimate = y_mean, std.error = sigma / sqrt(n)),
      est = tibble(.fitted = fits, .resid = res),
      fit = tibble(sigma2 = sigma^2),
      spec = tibble()
    ),
    class = "model_mean"
  )
}

#' Mean models
#' 
#' \code{MEAN()} returns an iid model applied to the formula's response variable.
#' 
#' The model does not use any specials, and so everything on the formula's 
#' right-hand-side will be ignored.
#' 
#' @param formula Model specification.
#' @param ... Not used.
#' 
#' @section Specials:
#' 
#' This model does not support usage of any specials. It only computes the mean!
#' 
#' @examples 
#' library(tsibbledata)
#' aus_elec %>% 
#'   model(rw = MEAN(Demand))
#' 
#' @export
MEAN <- function(formula, ...){
  mean_model <- new_model_class("mean", train = train_mean, specials = NULL)
  new_model_definition(mean_model, !!enquo(formula), ...)
}

#' @importFrom fablelite forecast
#' @importFrom stats qnorm time
#' @importFrom utils tail
#' @export
forecast.model_mean <- function(object, new_data, specials = NULL, bootstrap = FALSE, times = 5000, ...){
  h <- NROW(new_data)
  
  y_mean <- object$par$estimate
  n <- NROW(object$est)
  sigma <- sqrt(object$fit$sigma2)
  
  # Point forecasts
  fc <- rep(y_mean, h)
  
  # Intervals
  if (bootstrap){ # Compute prediction intervals using simulations
    sim <- map(seq_len(times), function(x){
      imitate(object, new_data, bootstrap = TRUE)[[".sim"]]
    }) %>%
      transpose %>%
      map(as.numeric)
    se <- map_dbl(sim, stats::sd)
    dist <- dist_sim(sim)
  }  else {
    se <- sigma * sqrt(1 + 1 / n)
    dist <- dist_normal(fc, se)
  }
  
  construct_fc(fc, se, dist)
}

#' @importFrom stats na.omit
#' @export
imitate.model_mean <- function(object, new_data, bootstrap = FALSE, ...){
  res <- residuals(object)
  f <- object$par$estimate
  
  if(is.null(new_data[[".innov"]])){
    if(bootstrap){
      new_data[[".innov"]] <- sample(na.omit(res) - mean(res, na.rm = TRUE),
                                     NROW(new_data), replace = TRUE)
    }
    else{
      new_data[[".innov"]] <- stats::rnorm(NROW(new_data), 
                                           sd = sqrt(object$fit$sigma2))
    }
  }
  
  new_data %>% 
    group_by_key() %>% 
    transmute(".sim" := f + !!sym(".innov"))
}


#' @export
fitted.model_mean <- function(object, ...){
  object$est[[".fitted"]]
}

#' @export
residuals.model_mean <- function(object, ...){
  object$est[[".resid"]]
}

#' @export
glance.model_mean <- function(x, ...){
  x$fit
}

#' @export
tidy.model_mean <- function(x, ...){
  x$par
}

#' @export
report.model_mean <- function(object, ...){
  cat("\n")
  cat(paste("Mean:", round(object$par$estimate, 4), "\n"))
  cat(paste("sigma^2:", round(object$fit$sigma2, 4), "\n"))
}

#' @export
model_sum.model_mean <- function(x){
  paste0("MEAN")#, ", intToUtf8(0x3BC), "=", format(x$par$estimate))
}