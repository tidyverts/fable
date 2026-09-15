context("test-lm.R")

test_that("LM", {
  skip_if_not_installed("forecast")

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

test_that("stream.TSLM", {
  library(tsibble)
  lung_deaths_male <- as_tsibble(mdeaths)

  train <- lung_deaths_male %>% filter(index < yearmonth("1979 Jan"))
  new_obs <- lung_deaths_male %>% filter(index >= yearmonth("1979 Jan"))

  # Streaming should exactly match a fixed-coefficient refit
  expect_stream_matches_refit <- function(streamed, refitted, fit) {
    # Coefficients are not re-estimated by stream()
    expect_equal(tidy(streamed)$estimate, tidy(fit)$estimate)
    expect_equal(tidy(streamed), tidy(refitted))
    expect_equal(fitted(streamed), fitted(refitted))
    expect_equal(residuals(streamed), residuals(refitted))
    expect_equal(glance(streamed), glance(refitted))
  }

  # trend() + season()
  fit_ts <- train %>% model(lm = TSLM(value ~ trend() + season()))
  streamed_ts <- fit_ts %>% stream(new_obs)
  refit_ts <- fit_ts %>% refit(lung_deaths_male, reestimate = FALSE)
  expect_stream_matches_refit(streamed_ts, refit_ts, fit_ts)
  expect_equal(NROW(fitted(streamed_ts)), NROW(lung_deaths_male))

  # External xreg
  train_x <- UKLungDeaths %>% filter(index < yearmonth("1979 Jan"))
  new_obs_x <- UKLungDeaths %>% filter(index >= yearmonth("1979 Jan"))
  fit_x <- train_x %>% model(lm = TSLM(mdeaths ~ trend() + fdeaths))
  streamed_x <- fit_x %>% stream(new_obs_x)
  refit_x <- fit_x %>% refit(UKLungDeaths, reestimate = FALSE)
  expect_stream_matches_refit(streamed_x, refit_x, fit_x)

  # Rank-deficient xreg
  UKLungDeaths_rd <- UKLungDeaths
  UKLungDeaths_rd$fdeaths2 <- UKLungDeaths_rd$fdeaths * 2
  fit_rd <- UKLungDeaths_rd %>%
    filter(index < yearmonth("1979 Jan")) %>%
    model(lm = TSLM(mdeaths ~ fdeaths + fdeaths2))
  streamed_rd <- fit_rd %>%
    stream(UKLungDeaths_rd %>% filter(index >= yearmonth("1979 Jan")))
  refit_rd <- fit_rd %>% refit(UKLungDeaths_rd, reestimate = FALSE)
  expect_stream_matches_refit(streamed_rd, refit_rd, fit_rd)

  # Chained single-observation stream() calls should match one large stream()
  chained_ts <- fit_ts
  for (i in seq_len(NROW(new_obs))) {
    chained_ts <- chained_ts %>% stream(new_obs[i, ])
  }
  expect_stream_matches_refit(chained_ts, streamed_ts, fit_ts)

  # Forecasts (including intervals) should match the refitted model
  fc_streamed <- streamed_ts %>% forecast(h = 12)
  fc_refit <- refit_ts %>% forecast(h = 12)
  expect_equal(fc_streamed$value, fc_refit$value)
  expect_equal(
    forecast(streamed_ts, h = 12, approx_normal = FALSE)$value,
    forecast(refit_ts, h = 12, approx_normal = FALSE)$value
  )

  # Streaming must start immediately after the trained data
  expect_error(
    fit_ts %>% stream(lung_deaths_male %>% filter(index >= yearmonth("1979 Feb"))),
    "must start one step beyond the end of"
  )
})
