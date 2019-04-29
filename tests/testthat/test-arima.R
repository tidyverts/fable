context("test-arima.R")

test_that("ARIMA", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% model(arima = ARIMA(value))
  stats_fit <- arima(USAccDeaths, c(0,1,1), list(order = c(0,1,1), 12))
  
  expect_identical(
    coef(fable_fit$arima[[1]]$fit$model),
    coef(stats_fit)
  )
  
  # Automatic d/D selection
  fit <- USAccDeaths_tbl %>% model(ARIMA(value ~ pdq(q=1)))
  expect_identical(
    model_sum(fit[[1]][[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
  
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(model = ARIMA(value ~ pdq(0,1,1) + PDQ(0,1,1)))

  expect_identical(
    coef(fable_fit$model[[1]]$fit$model),
    coef(stats_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
  
  fable_fc <- fable_fit %>% forecast
  stats_fc <- stats_fit %>% predict(24)
  
  expect_equivalent(
    fable_fc$value,
    unclass(stats_fc$pred)
  )
})


test_that("ARIMA with xregs", {
  tr <- UKLungDeaths %>% head(-12)
  ts <- UKLungDeaths %>% tail(12)
  fable_fit <- tr %>% model(model = ARIMA(mdeaths ~ fdeaths + pdq(d=1) + PDQ(0,0,0)))
  stats_fit <- arima(head(mdeaths, -12), c(0,1,1),
                     xreg = data.frame(fdeaths = head(fdeaths, -12)))
  
  expect_equivalent(
    coef(fable_fit$model[[1]]$fit$model),
    coef(stats_fit)
  )
  
  fable_fc <- fable_fit %>% forecast(ts)
  stats_fc <- stats_fit %>% predict(12, newxreg = data.frame(fdeaths = tail(fdeaths, 12)))
  
  expect_equivalent(
    fable_fc$mdeaths,
    unclass(stats_fc$pred)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "LM w/ ARIMA(0,1,1) errors"
  )
})