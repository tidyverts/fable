context("test-arima.R")

test_that("ARIMA", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% ARIMA(value ~ pdq(d = 1) + PDQ(D = 1))
  stats_fit <- arima(USAccDeaths, c(0,1,1), list(order = c(0,1,1), 12))
  
  expect_identical(
    coef(fable_fit$model[[1]]$model),
    coef(stats_fit)
  )
  
  # Lack of support for automatic differencing
  expect_error(
    USAccDeaths_tbl %>% ARIMA(value ~ pdq(q=1)),
    "Automatic selection of differencing is currently not implemented"
  )
  
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% ARIMA(value ~ pdq(0,1,1) + PDQ(0,1,1))
  
  expect_identical(
    coef(fable_fit$model[[1]]$model),
    coef(stats_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
  
  fable_fc <- fable_fit %>% forecast
  stats_fc <- stats_fit %>% predict(24)
  
  expect_equivalent(
    fable_fc$mean,
    unclass(stats_fc$pred)
  )
})


test_that("ARIMA with xregs", {
  tr <- UKLungDeaths %>% head(-12)
  ts <- UKLungDeaths %>% tail(12)
  fable_fit <- tr %>% ARIMA(mdeaths ~ fdeaths + pdq(d=1) + PDQ(D=1))
  forecast_fit <- forecast::auto.arima(head(mdeaths, -12), xreg = head(fdeaths, -12))
  
  expect_equivalent(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  fable_fc <- fable_fit %>% forecast(12, xreg = ts[["fdeaths"]])
  forecast_fc <- forecast_fit %>% forecast(xreg = tail(fdeaths, 12))
  
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ARIMA(0,1,1)"
  )
})