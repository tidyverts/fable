context("test-ets.R")

forecast_fit <- USAccDeaths %>% forecast::ets()
test_that("Automatic ETS selection", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value))
  
  expect_equivalent(
    tidy(fable_fit$ets[[1]]$fit)$estimate,
    coef(forecast_fit)
  )
  
  # Short series
  expect_equal(
    tidy(model(UKLungDeaths[1:24,], ETS(mdeaths)))$estimate,
    c(1, 2134),
    tolerance = 0.5
  )
})

test_that("Manual ETS selection", {
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value ~ error("A") + trend("N") + season("A")))
  
  expect_equivalent(
    tidy(fable_fit$ets[[1]]$fit)$estimate,
    coef(forecast_fit)
  )
  
  expect_identical(
    model_sum(fable_fit$ets[[1]]),
    "ETS(A,N,A)"
  )
  
  fable_fc <- fable_fit %>% forecast
  forecast_fc <- forecast_fit %>% forecast::forecast()
  
  expect_equivalent(
    fable_fc$value,
    unclass(forecast_fc$mean)
  )
  
  # Test simulation
  fable_fit %>%
    imitate(USAccDeaths_tbl)
  fable_fit %>% 
    imitate(USAccDeaths_tbl %>% 
              dplyr::mutate(index = index + 72))
  
  # Test refit
  expect_identical(
    tidy(refit(fable_fit, USAccDeaths_tbl))$estimate == tidy(fable_fit)$estimate,
    c(rep(TRUE, 2), rep(FALSE, 12))
  )
  expect_identical(
    tidy(refit(fable_fit, USAccDeaths_tbl, reinitialise = FALSE))$estimate,
    tidy(fable_fit)$estimate
  )
  
  # Test components
  cmp <- components(fable_fit)
  expect_identical(
    tidy(fable_fit)$estimate[3:14],
    c(cmp$level[12], cmp$season[12:2])
  )
  expect_s3_class(
    cmp, "dcmp_ts"
  )
  
  # Test report
  expect_output(
    report(fable_fit),
    "sigma:  292.6907"
  )
  
  aug <- augment(fable_fit)
  expect_equal(
    aug$value,
    aug$.fitted + aug$.resid
  )
})


test_that("ETS with bad inputs", {
  # Test for multiple specials defined
  expect_error(
    USAccDeaths_tbl %>% model(ETS(value ~ error("A") + error("A"))),
    "Only one special of each type is allowed for ETS"
  )
  
  expect_error(
    UKLungDeaths %>% 
      model(ETS(vars(mdeaths, fdeaths))),
    "Only univariate responses are supported by ETS"
  )
  
  UK_missing <- UKLungDeaths
  UK_missing[["mdeaths"]][3:5] <- NA
  expect_error(
    UK_missing %>% 
      model(ETS(mdeaths)),
    "ETS does not support missing values"
  )
  
  expect_error(
    UKLungDeaths %>% 
      model(ETS(mdeaths ~ trend("M") + season("A"))),
    "No valid ETS models have been allowed"
  )
  
  expect_error(
    UKLungDeaths[1:2,] %>% 
      model(ETS(mdeaths)),
    "Not enough data to estimate this ETS model"
  )
})


test_that("Multiplicative ETS models", {
  fable_fit <- USAccDeaths_tbl %>%
    model(ets = ETS(value ~ error("M") + trend("N") + season("N")))
  expect_true(
    is.constant(forecast(fable_fit)$value)
  )
  
  expect_s3_class(
    USAccDeaths_tbl %>%
      model(ets = ETS(value ~ error("M") + trend("A") + season("M"))) %>% 
      forecast(),
    "fbl_ts"
  )
  
  
  expect_s3_class(
    USAccDeaths_tbl %>%
      model(ets = ETS(value ~ error("M") + trend("M") + season("M"))) %>% 
      forecast(times = 5),
    "fbl_ts"
  )
})