context("test-lm.R")

test_that("LM", {
  # NULL model selection
  fable_fit <- USAccDeaths_tbl %>% LM(value)
  forecast_fit <- lm(USAccDeaths~1)
  
  expect_identical(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  # Trend + Season
  fable_fit <- USAccDeaths_tbl %>% LM(value ~ trend() + season())
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + season)
  
  expect_equivalent(
    unclass(fitted(fable_fit$model[[1]])),
    unclass(fitted(forecast_fit))
  )
  
  # Fourier
  fable_fit <- USAccDeaths_tbl %>% LM(value ~ trend() + fourier(K=5))
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + forecast::fourier(USAccDeaths, K=5))
  
  expect_equivalent(
    unclass(fitted(fable_fit$model[[1]])),
    unclass(fitted(forecast_fit))
  )
  
  # Model summary
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "LM"
  )
  
  # Forecast
  fable_fc <- fable_fit %>% forecast()
  forecast_fc <- forecast_fit %>% forecast(h=12)
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
})
