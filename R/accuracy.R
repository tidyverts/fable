ME <- function(res, na.rm = TRUE, ...){
  mean(res, na.rm = na.rm, ...)
}

MSE <- function(res, na.rm = TRUE, ...){
  mean(res ^ 2, na.rm = na.rm, ...)
}

RMSE <- function(res, na.rm = TRUE, ...){
  sqrt(MSE(res, na.rm = na.rm, ...))
}

MAE <- function(res, na.rm = TRUE, ...){
  mean(abs(res), na.rm = na.rm, ...)
}

MPE <- function(res, y, na.rm = TRUE, ...){
  mean(res / y * 100, na.rm = TRUE, ...)
}

MAPE <- function(res, y, na.rm = TRUE, ...){
  mean(abs(res / y * 100), na.rm = TRUE, ...)
}

MASE <- function(res, x, demean = FALSE, na.rm = TRUE, ...){
  if(demean){
    scale <- mean(abs(x - mean(x, na.rm = na.rm, ...)), na.rm = na.rm, ...)
  }
  else{
    scale <- mean(abs(x), na.rm = na.rm, ...)
  }
  mase <- mean(abs(res / scale), na.rm = na.rm, ...)
}

ACF1 <- function(res, na.action = stats::na.pass, ...){
  stats::acf(res, plot = FALSE, lag.max = 2, na.action = na.action, ...)$acf[2, 1, 1]
}