context("test-ets.R")
skip_if_not_installed("forecast")

forecast_fit <- USAccDeaths %>% forecast::ets()
test_that("Automatic ETS selection", {
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value))

  expect_equivalent(
    tidy(fable_fit$ets[[1]]$fit)$estimate,
    c(coef(forecast_fit), -sum(coef(forecast_fit)[-(1:3)]))
  )

  # Short series
  expect_equal(
    tidy(model(UKLungDeaths[1:24, ], ETS(mdeaths)))$estimate,
    c(1, 2134),
    tolerance = 0.5
  )
})

test_that("Manual ETS selection", {
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(ets = ETS(value ~ error("A") + trend("N") + season("A")))

  expect_equivalent(
    tidy(fable_fit$ets[[1]]$fit)$estimate,
    c(coef(forecast_fit), -sum(coef(forecast_fit)[-(1:3)]))
  )

  expect_identical(
    model_sum(fable_fit$ets[[1]]),
    "ETS(A,N,A)"
  )

  fable_fc <- fable_fit %>% forecast()
  forecast_fc <- forecast_fit %>% forecast::forecast()

  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(forecast_fc$mean)
  )

  # Test simulation
  fable_fit %>%
    generate(USAccDeaths_tbl)
  fable_fit %>%
    generate(USAccDeaths_tbl %>%
      dplyr::mutate(index = index + 72))

  # Test refit
  expect_identical(
    tidy(refit(fable_fit, USAccDeaths_tbl))$estimate == tidy(fable_fit)$estimate,
    c(rep(TRUE, 2), rep(FALSE, 13))
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
    "sigma\\^2:  85667.86"
  )

  aug <- augment(fable_fit)
  expect_equal(
    aug$value,
    aug$.fitted + aug$.resid
  )

  # Test specification of smoothing params
  coef <- USAccDeaths_tbl %>%
    model(ETS(value ~ error("A") + season("A", gamma = 0.0001) +
      trend("Ad", alpha = 0.5, beta = 0.006, phi = 0.975))) %>%
    tidy()
  expect_identical(
    coef$estimate[1:4],
    c(0.5, 0.006, 0.0001, 0.975)
  )
  expect_identical(
    coef$term,
    c("alpha", "beta", "gamma", "phi", "l[0]", "b[0]", sprintf("s[%i]", 0:-11))
  )
})


test_that("ETS with bad inputs", {
  # Test for multiple specials defined
  expect_warning(
    USAccDeaths_tbl %>% model(ETS(value ~ error("A") + error("A"))),
    "Only one special of each type is allowed for ETS"
  )

  expect_warning(
    USAccDeaths_tbl %>% model(ETS(value ~ trend(alpha = 1.5))),
    "Inconsistent parameter boundaries"
  )

  expect_warning(
    USAccDeaths_tbl %>% model(ETS(value ~ error("A") + trend("A", alpha = 0.2, beta = 0.5) + season("N"))),
    "Parameters out of range"
  )

  expect_warning(
    UKLungDeaths %>%
      model(ETS(vars(mdeaths, fdeaths))),
    "Only univariate responses are supported by ETS"
  )

  UK_missing <- UKLungDeaths
  UK_missing[["mdeaths"]][3:5] <- NA
  expect_warning(
    UK_missing %>%
      model(ETS(mdeaths)),
    "ETS does not support missing values"
  )

  expect_warning(
    UKLungDeaths %>%
      model(ETS(mdeaths ~ trend("M") + season("A"))),
    "No valid ETS models have been allowed"
  )

  expect_warning(
    UKLungDeaths[1:2, ] %>%
      model(ETS(mdeaths)),
    "Not enough data to estimate this ETS model"
  )
})


test_that("Multiplicative ETS models", {
  fable_fit <- USAccDeaths_tbl %>%
    model(ets = ETS(value ~ error("M") + trend("N") + season("N")))
  expect_true(
    is.constant(fc_mean(forecast(fable_fit)$value))
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

test_that("Automatic ETS selection bug (#425)", {
  train <- tsibble(
    YM = yearmonth("2022 Jan") + 0:35,
    value = rep(c(rep(-78040, 11), -78061), 3),
    index = YM
  )
  expect_identical(
    model_sum(model(train, ets=ETS(value))$ets[[1]]),
    "ETS(A,N,N)"
  )
})