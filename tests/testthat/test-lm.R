context("test-lm.R")

test_that("LM", {
  # NULL model selection
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value))
  forecast_fit <- lm(USAccDeaths ~ 1)

  expect_equivalent(
    coef(fable_fit)$estimate,
    coef(forecast_fit)
  )

  # Trend + Season
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value ~ trend() + season()))
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + season)

  expect_equivalent(
    unclass(fitted(fable_fit)[[".fitted"]]),
    unclass(fitted(forecast_fit))
  )

  # Model coefs
  expect_equivalent(
    tidy(fable_fit) %>% dplyr::filter(term == "trend()") %>% dplyr::pull(estimate),
    coef(forecast_fit)["trend"]
  )

  # Forecast
  fable_fc <- fable_fit %>% forecast(h = 12)
  fable_fc_short <- fable_fit %>% forecast(h = 1)
  forecast_fc <- forecast_fit %>% forecast::forecast(h = 12)
  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(forecast_fc$mean)
  )
  expect_equivalent(
    fable_fc$value[1],
    fable_fc_short$value
  )

  fable_fc_sim <- fable_fit %>% forecast(h = 12, bootstrap = TRUE, times = 5)
  # expect_equal(
  #   fc_mean(fable_fc$value),
  #   fc_mean(fable_fc_sim$value)
  # )

  # Fourier
  fable_fit <- USAccDeaths_tbl %>% model(lm = TSLM(value ~ trend() + fourier(K = 5)))
  forecast_fit <- forecast::tslm(USAccDeaths ~ trend + forecast::fourier(USAccDeaths, K = 5))

  expect_equivalent(
    unclass(fitted(fable_fit)[[".fitted"]]),
    unclass(fitted(forecast_fit))
  )

  # Model summary
  expect_identical(
    model_sum(fable_fit$lm[[1]]),
    "TSLM"
  )

  # Model report
  expect_output(
    report(fable_fit),
    "Residual standard error: 442.5"
  )

  # Model glance
  expect_equal(
    with(glance(fable_fit), df + df.residual),
    NROW(USAccDeaths_tbl)
  )

  # Refit
  expect_identical(
    tidy(fable_fit)$estimate,
    tidy(refit(fable_fit, USAccDeaths_tbl))$estimate
  )

  # Interpolate
  USAccDeaths_tbl[["value"]][10] <- NA
  expect_equal(
    interpolate(fable_fit, USAccDeaths_tbl)[["value"]][10],
    fitted(fable_fit)[[".fitted"]][10]
  )
})
