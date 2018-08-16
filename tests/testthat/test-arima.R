context("test-arima.R")

test_that("ARIMA", {
  # Automatic model selection
  fit <- USAccDeaths_tbl %>% ARIMA(value)
  fc_fit <- USAccDeaths %>% forecast::auto.arima()
  
  expect_identical(
    fit$model[[1]]$coef,
    fc_fit$coef
  )
  
  # Partial automatic model selection
  expect_message(
    fit <- USAccDeaths_tbl %>% ARIMA(value ~ pdq(q=1)),
    "Partial automation of parameters is not yet supported"
  )
  
  # Manual model selection
  fit <- USAccDeaths_tbl %>% ARIMA(value ~ pdq(0,1,1) + PDQ(0,1,1))
  
  expect_identical(
    fit$model[[1]]$coef,
    fc_fit$coef
  )
  
  expect_identical(
    model_sum(fit$model[[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
})