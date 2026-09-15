context("test-rw.R")
test_that("NAIVE", {
  skip_if_not_installed("forecast")
  fable_fit <- USAccDeaths_tbl %>% model(naive = NAIVE(value))
  forecast_fc <- forecast::naive(USAccDeaths, h = 12)

  expect_equivalent(
    fitted(fable_fit)[[".fitted"]],
    unclass(fitted(forecast_fc))
  )

  fable_fc <- fable_fit %>% forecast(h = 12)

  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(forecast_fc$mean)
  )
  expect_identical(
    model_sum(fable_fit$naive[[1]]),
    "NAIVE"
  )
})

test_that("RW w/ drift", {
  skip_if_not_installed("forecast")
  fable_fit <- USAccDeaths_tbl %>% model(rw = RW(value ~ drift()))
  forecast_fc <- forecast::rwf(USAccDeaths, drift = TRUE, h = 12)

  expect_equivalent(
    fitted(fable_fit)[[".fitted"]],
    unclass(fitted(forecast_fc))
  )

  fable_fc <- fable_fit %>% forecast(h = 12)

  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(forecast_fc$mean)
  )

  if(packageVersion("forecast") > "8.17.0"){
    expect_equivalent(
      unclass(hilo(fable_fc)$`80%`)$upper,
      unclass(forecast_fc$upper[, 1])
    )
  }

  expect_identical(
    model_sum(fable_fit$rw[[1]]),
    "RW w/ drift"
  )

  expect_output(
    report(fable_fit),
    "Drift: 3\\.2817 \\(se: 87\\.2696\\)"
  )

  expect_equal(
    tidy(fable_fit)$estimate,
    forecast_fc$model$par$drift
  )

  expect_equal(
    glance(fable_fit)$sigma2,
    forecast_fc$model$sigma2
  )

  expect_equivalent(
    residuals(fable_fit)[[".resid"]],
    unclass(residuals(forecast_fc))
  )
})

test_that("SNAIVE", {
  skip_if_not_installed("forecast")
  fable_fit <- USAccDeaths_tbl %>% model(snaive = SNAIVE(value))
  forecast_fc <- forecast::snaive(USAccDeaths, h = 12)

  expect_equivalent(
    fitted(fable_fit)[[".fitted"]],
    unclass(fitted(forecast_fc))
  )

  fable_fc <- fable_fit %>% forecast(h = 12)

  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(forecast_fc$mean)
  )

  expect_identical(
    model_sum(fable_fit$snaive[[1]]),
    "SNAIVE"
  )

  fable_fc_sim <- fable_fit %>%
    forecast(h = 12, bootstrap = TRUE, times = 5)
  # expect_equal(
  #   fable_fc$value,
  #   fable_fc_sim$value
  # )

  skip_if_not_installed("tsibbledata")
  expect_warning(
    tsibbledata::vic_elec %>%
      tsibble::index_by(date = as.Date(Time)) %>%
      dplyr::summarise(demand = mean(Demand)) %>%
      model(SNAIVE(demand ~ lag("year"))),
    "Non-integer lag orders for random walk models are not supported"
  )
})

test_that("RW short", {
  library(tsibble)
  fc <- suppressWarnings(tsibble(year = 2010:2012, y = 1:3, index = year) %>%
    model(SNAIVE(y ~ lag(4))) %>%
    forecast(h = 4))

  expect_equal(
    fc_mean(fc$y),
    c(NA, 1, 2, 3)
  )
})

test_that("stream.RW", {
  library(tsibble)
  lung_deaths_male <- as_tsibble(mdeaths)

  train <- lung_deaths_male %>% filter(index < yearmonth("1979 Jan"))
  new_obs <- lung_deaths_male %>% filter(index >= yearmonth("1979 Jan"))

  # SNAIVE: streaming should exactly match a fixed-coefficient refit
  fit_sn <- train %>% model(snaive = SNAIVE(value))
  streamed_sn <- fit_sn %>% stream(new_obs)
  refit_sn <- fit_sn %>% refit(lung_deaths_male, reestimate = FALSE)

  expect_equal(
    fitted(streamed_sn)[[".fitted"]],
    fitted(refit_sn)[[".fitted"]]
  )
  expect_equal(
    residuals(streamed_sn)[[".resid"]],
    residuals(refit_sn)[[".resid"]]
  )
  expect_equal(glance(streamed_sn)$sigma2, glance(refit_sn)$sigma2)

  # RW w/ drift: the drift coefficient should not be re-estimated by stream()
  fit_rw <- train %>% model(rw = RW(value ~ drift()))
  streamed_rw <- fit_rw %>% stream(new_obs)
  refit_rw <- fit_rw %>% refit(lung_deaths_male, reestimate = FALSE)

  expect_equal(tidy(streamed_rw)$estimate, tidy(fit_rw)$estimate)
  expect_equal(
    fitted(streamed_rw)[[".fitted"]],
    fitted(refit_rw)[[".fitted"]]
  )
  expect_equal(
    residuals(streamed_rw)[[".resid"]],
    residuals(refit_rw)[[".resid"]]
  )
  expect_equal(glance(streamed_rw)$sigma2, glance(refit_rw)$sigma2)

  # Chained stream() calls should give the same result as one large stream()
  mid <- lung_deaths_male %>% filter(index < yearmonth("1978 Jul"))
  part1 <- lung_deaths_male %>%
    filter(index >= yearmonth("1978 Jul"), index < yearmonth("1979 Jan"))
  part2 <- new_obs

  chained_sn <- mid %>%
    model(snaive = SNAIVE(value)) %>%
    stream(part1) %>%
    stream(part2)

  expect_equal(
    fitted(chained_sn)[[".fitted"]],
    fitted(refit_sn)[[".fitted"]]
  )

  # Streaming must start immediately after the trained data
  expect_error(
    fit_sn %>% stream(lung_deaths_male %>% filter(index >= yearmonth("1979 Feb"))),
    "must start one step beyond the end of"
  )
})

test_that("lagwalk with bad inputs", {
  expect_warning(
    UKLungDeaths %>%
      model(SNAIVE(vars(mdeaths, fdeaths))),
    "Only univariate responses are supported by lagwalks"
  )

  expect_warning(
    UKLungDeaths %>%
      model(SNAIVE(resp(rlang::rep_along(mdeaths, NA)))),
    "All observations are missing"
  )

  expect_warning(
    UKLungDeaths %>%
      model(SNAIVE(mdeaths ~ lag(1))),
    "Non-seasonal model specification provided"
  )
})
