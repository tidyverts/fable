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

MASE <- function(res, y, demean = FALSE, na.rm = TRUE, period, d = period == 1, D = period > 1, ...){
  if (D > 0) { # seasonal differencing
    y <- diff(y, lag = period, differences = D)
  }
  if (d > 0) {
    y <- diff(y, differences = d)
  }
  if(demean){
    scale <- mean(abs(y - mean(y, na.rm = na.rm, ...)), na.rm = na.rm, ...)
  }
  else{
    scale <- mean(abs(y), na.rm = na.rm, ...)
  }
  mase <- mean(abs(res / scale), na.rm = na.rm, ...)
}

ACF1 <- function(res, na.action = stats::na.pass, ...){
  stats::acf(res, plot = FALSE, lag.max = 2, na.action = na.action, ...)$acf[2, 1, 1]
}

accuracy.mable <- function(f, period = "smallest"){
  accuracy_data <- residuals(f) %>% left_join(getResponse(f), by = expr_text(index(.)))
  period <- get_frequencies(period, accuracy_data)
  accuracy_data %>% 
    group_by(!!!syms(key_vars(.))) %>% 
    as_tibble() %>% 
    summarise(
      Type = "Training set",
      ME = ME(residuals),
      RMSE = RMSE(residuals),
      MAE = MAE(residuals),
      MPE = MPE(residuals, response),
      MAPE = MAPE(residuals, response),
      MASE = MASE(residuals, response, period = period),
      ACF1 = ACF1(residuals)
    )
}

accuracy.fable <- function(f, x, period = "smallest"){
  response <- getResponse(f)$response
  browser()
  accuracy_data <- residuals(f) %>% mutate(Type = "Training set")
  if(!missing(x)){
    keys <- key_vars(f)
    
    test_data <- x %>% 
      group_by(!!!syms(keys)) %>% 
      nest(.key = "newdata")
    
    if(NROW(f) > 1){
      left_join(f, test_data, by = keys)
    }
    else{
      f$newdata <- test_data$newdata
    }
    
    resid_test <- f %>% 
      transmute(!!!syms(keys),
                residuals = pmap(list(model, forecast, newdata),
                                 function(model, forecast, newdata){
                                   forecast %>% 
                                     transmute(residuals = newdata[[model%@%"response"]] - mean)
                                 })
      ) %>% 
      unnest %>% 
      mutate(Type = "Test set")
    
    accuracy_data <- rbind(accuracy_data, resid_test)
  }
  period <- get_frequencies(period, accuracy_data)
  accuracy_data %>% 
    group_by(!!!syms(c(key_vars(.), "Type"))) %>% 
    as_tibble() %>% 
    summarise(
      ME = ME(residuals),
      RMSE = RMSE(residuals),
      MAE = MAE(residuals),
      MPE = MPE(residuals, response),
      MAPE = MAPE(residuals, response),
      MASE = MASE(residuals, response, period = period),
      ACF1 = ACF1(residuals)
    )
}