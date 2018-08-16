context("test-arima.R")

test_that("ARIMA", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% ARIMA(value)
  forecast_fit <- USAccDeaths %>% forecast::auto.arima()
  
  expect_identical(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  # Partial automatic model selection
  expect_message(
    USAccDeaths_tbl %>% ARIMA(value ~ pdq(q=1)),
    "Partial automation of parameters is not yet supported"
  )
  
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% ARIMA(value ~ pdq(0,1,1) + PDQ(0,1,1))
  
  expect_identical(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
  
  fable_fc <- fable_fit %>% forecast
  forecast_fc <- forecast_fit %>% forecast
  
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
})
