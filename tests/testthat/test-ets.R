context("test-ets.R")

test_that("ETS", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value))
  forecast_fit <- USAccDeaths %>% forecast::ets()
  
  expect_identical(
    coef(fable_fit$ets[[1]]$fit),
    coef(forecast_fit)
  )
  
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value ~ error("A") + trend("N") + season("A")))
  
  expect_identical(
    coef(fable_fit$ets[[1]]$fit),
    coef(forecast_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$ets[[1]]),
    "ETS(A,N,A)"
  )
  
  # Test for multiple specials defined
  expect_error(
    USAccDeaths_tbl %>% model(ETS(value ~ error("A") + error("A"))),
    "Only one special of each type is allowed for ETS."
  )
  
  fable_fc <- fable_fit %>% forecast
  forecast_fc <- forecast_fit %>% forecast
  
  expect_equivalent(
    fable_fc$value,
    unclass(forecast_fc$mean)
  )
  
  # Test simulation
  fable_fit %>% 
    simulate(USAccDeaths_tbl)
  fable_fit %>% 
    simulate(USAccDeaths_tbl %>% 
               tsibble::mutate(index = index + 72))
})
