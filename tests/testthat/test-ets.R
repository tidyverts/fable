context("test-ets.R")

test_that("ETS", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% ETS(value)
  forecast_fit <- USAccDeaths %>% forecast::ets()
  
  expect_identical(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% ETS(value ~ error("A") + trend("N") + season("A"))
  
  expect_identical(
    coef(fable_fit$model[[1]]),
    coef(forecast_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ETS(A,N,A)"
  )
  
  # Test for multiple specials defined
  expect_error(
    USAccDeaths_tbl %>% ETS(value ~ error("A") + error("A")),
    "Only one special of each type is allowed for ETS."
  )
})
