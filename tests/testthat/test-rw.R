context("test-rw.R")
test_that("NAIVE", {
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

  if(packageVersion("fabletools") > "0.1.3"){
    expect_equivalent(
      unclass(hilo(fable_fc)$`80%`)$upper,
      unclass(forecast_fc$upper[, 1])
    )
  } else {
    expect_equivalent(
      hilo(fable_fc)$`80%`$.upper,
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
