skip_if_not_installed("mixtime")
skip_if_not(
  inherits(tsibble::interval(tsibble::tsibble(t = 1:2, y = 1, index = t)), "mixtime::mt_unit"),
  "tsibble mixtime intervals not available"
)

set.seed(1)
y <- rnorm(72) + 10
mt <- tsibble::tsibble(t = mixtime::yearmonth(540:611), y = y, index = t)
legacy <- tsibble::tsibble(t = tsibble::yearmonth(540:611), y = y, index = t)

test_that("ETS supports monthly mixtime indices", {
  fit <- model(mt, ets = ETS(y))
  fit_legacy <- model(legacy, ets = ETS(y))

  expect_equal(forecast(fit, h = 6)$.mean, forecast(fit_legacy, h = 6)$.mean)
  expect_equal(NROW(generate(fit, h = 6)), 6)
  expect_s3_class(refit(fit, mt)$ets[[1]], "mdl_ts")
  expect_equal(NROW(components(fit)), NROW(components(fit_legacy)))
})

test_that("ETS supports hourly mixtime seasonal periods", {
  hourly <- tsibble::tsibble(
    t = mixtime::datetime(seq(as.POSIXct("2020-01-01", tz = "UTC"), by = "hour", length.out = 96)),
    y = rnorm(96) + 10, index = t
  )
  fit <- model(hourly, ets = ETS(y ~ season(period = "day")))
  expect_equal(NROW(forecast(fit, h = 12)), 12)
})

test_that("RW and SNAIVE bootstrap generation supports mixtime indices", {
  fit <- model(mt, rw = RW(y), snaive = SNAIVE(y))
  expect_equal(NROW(generate(fit, h = 6, bootstrap = TRUE)), 12)
})

test_that("stream() checks the start of mixtime new_data", {
  fit <- model(mt, rw = RW(y), arima = ARIMA(y), mean = MEAN(y))
  bad <- tsibble::tsibble(t = mixtime::yearmonth(615:616), y = 10, index = t)
  expect_error(stream(fit["rw"], bad), "one step beyond")
  expect_error(stream(fit["arima"], bad), "one step beyond")
  expect_error(stream(fit["mean"], bad), "one step beyond")
})
