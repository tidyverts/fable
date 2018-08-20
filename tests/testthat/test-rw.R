context("test-rw.R")

test_that("NAIVE", {
  fable_fit <- USAccDeaths_tbl %>% NAIVE(value)
  forecast_fc <- forecast::naive(USAccDeaths, h = 12)
  
  expect_equivalent(
    fitted(fable_fit)[["fitted"]],
    unclass(fitted(forecast_fc))
  )
  
  fable_fc <- fable_fit %>% forecast
  
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "NAIVE"
  )
})

test_that("RW w/ drift", {
  fable_fit <- USAccDeaths_tbl %>% RW(value ~ drift())
  forecast_fc <- forecast::rwf(USAccDeaths, drift = TRUE, h = 12)
  
  expect_equivalent(
    fitted(fable_fit)[["fitted"]],
    unclass(fitted(forecast_fc))
  )
  
  fable_fc <- fable_fit %>% forecast
  
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
  
  expect_equivalent(
    upper(summary(fable_fc)$`80%`),
    unclass(forecast_fc$upper[,1])
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "RW w/ drift"
  )
})

# test_that("RW w/ xreg", {
#   tr <- UKLungDeaths %>% head(-12)
#   ts <- UKLungDeaths %>% tail(12)
#   fable_fit <- tr %>% RW(mdeaths ~ fdeaths + drift())
#   
#   fable_fc <- fable_fit %>% forecast(newdata = ts)
#   
#   expect_identical(
#     model_sum(fable_fit$model[[1]]),
#     "lm w/ RW e"
#   )
# })

test_that("SNAIVE", {
  fable_fit <- USAccDeaths_tbl %>% SNAIVE(value)
  forecast_fc <- forecast::snaive(USAccDeaths, h = 12)
  
  expect_equivalent(
    fitted(fable_fit)[["fitted"]],
    unclass(fitted(forecast_fc))
  )
  
  fable_fc <- fable_fit %>% forecast
  
  expect_equivalent(
    summary(fable_fc)$mean,
    unclass(forecast_fc$mean)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "SNAIVE"
  )
})
