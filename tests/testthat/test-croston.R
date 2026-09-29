context("test-croston.R")

test_that("CROSTON forecasts from the final smoothed demand rate", {
  spec <- CROSTON(v ~ demand(initial = 3, param = 0.3) + interval(initial = 2, param = 0.3))
  v <- c(0, 3, 0, 0, 2, 0, 1, 0, 0, 4, 0, 5)
  ends_nonzero <- tsibble::tsibble(t = seq_along(v), v = v, index = t) %>% model(spec)
  ends_zero <- tsibble::tsibble(t = seq_len(13), v = c(v, 0), index = t) %>% model(spec)

  fc <- forecast(ends_nonzero, h = 2)$.mean
  expect_equal(fc, forecast(ends_zero, h = 2)$.mean)
  expect_equal(fc, rep(tail(fitted(ends_zero)$.fitted, 1), 2))
})
