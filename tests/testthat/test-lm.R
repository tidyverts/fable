context("test-lm.R")

test_that("LM", {
  # NULL model selection
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value))
  forecast_fit <- lm(USAccDeaths~1)
  
  expect_identical(
    coef(fable_fit$lm[[1]]$fit$model),
    coef(forecast_fit)
  )
  
  # Trend + Season
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value ~ trend() + season()))
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + season)
  
  expect_equivalent(
    unclass(fitted(fable_fit)[[".fitted"]]),
    unclass(fitted(forecast_fit))
  )
  
  # Forecast
  fable_fc <- fable_fit %>% forecast(h = 12)
  forecast_fc <- forecast_fit %>% forecast::forecast(h = 12)
  expect_equivalent(
    fable_fc$value,
    unclass(forecast_fc$mean)
  )
  
  # Fourier
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value ~ trend() + fourier(K=5)))
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + forecast::fourier(USAccDeaths, K=5))
  
  expect_equivalent(
    unclass(fitted(fable_fit)[[".fitted"]]),
    unclass(fitted(forecast_fit))
  )
  
  # Model summary
  expect_identical(
    model_sum(fable_fit$lm[[1]]),
    "TSLM"
  )
})
