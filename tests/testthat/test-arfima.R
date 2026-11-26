context("test-arfima.R")

sunspot_tbl <- as_tsibble(sunspot.year)

test_that("refit.fbl_ARFIMA keeps a single mu and d", {
  skip_if_not_installed("fracdiff")
  fit <- suppressWarnings(
    head(sunspot_tbl, -20) %>% model(arfima = ARFIMA(value ~ pdq(1, NULL, 1)))
  )
  refitted <- refit(fit, sunspot_tbl, reestimate = FALSE)

  expect_identical(tidy(refitted)$term, tidy(fit)$term)
  expect_identical(tidy(refitted)$estimate, tidy(fit)$estimate)
})

test_that("stream.fbl_ARFIMA matches a full refit", {
  skip_if_not_installed("fracdiff")
  tr <- head(sunspot_tbl, -20)
  nw <- tail(sunspot_tbl, 20)

  fit <- suppressWarnings(tr %>% model(arfima = ARFIMA(value ~ pdq(1, NULL, 1))))
  streamed <- fit %>% stream(nw)
  refitted <- fit %>% refit(sunspot_tbl, reestimate = FALSE)

  expect_s3_class(streamed$arfima[[1]]$fit, "fbl_ARFIMA")
  expect_identical(tidy(streamed)$estimate, tidy(fit)$estimate)
  expect_equal(NROW(augment(streamed)), NROW(sunspot_tbl))

  expect_equal(fitted(streamed)$.fitted, fitted(refitted)$.fitted)
  expect_equal(residuals(streamed)$.resid, residuals(refitted)$.resid)
  expect_equal(fitted(streamed)$.fitted[seq_len(NROW(tr))], fitted(fit)$.fitted)
  expect_equal(glance(streamed)$log_lik, glance(refitted)$log_lik)

  # sigma2 counts the trained ARMA coefficients, as in stream.ARIMA()
  expect_equal(
    glance(streamed)$sigma2,
    sum(residuals(streamed)$.resid^2) / (NROW(sunspot_tbl) - 3 + 1)
  )

  fc_streamed <- forecast(streamed, h = 10)
  fc_refitted <- forecast(refitted, h = 10)
  expect_equal(fc_streamed$.mean, fc_refitted$.mean)
  expect_equal(fc_streamed$value, fc_refitted$value, tolerance = 0.01)
})

test_that("Chained ARFIMA streaming matches one-shot streaming", {
  skip_if_not_installed("fracdiff")
  fit <- head(sunspot_tbl, -20) %>% model(arfima = ARFIMA(value ~ pdq(1, NULL, 0)))

  chained <- fit %>%
    stream(sunspot_tbl %>% tail(20) %>% head(10)) %>%
    stream(sunspot_tbl %>% tail(10))
  one_shot <- fit %>% stream(sunspot_tbl %>% tail(20))

  expect_equal(fitted(chained)$.fitted, fitted(one_shot)$.fitted)
  expect_equal(residuals(chained)$.resid, residuals(one_shot)$.resid)
  expect_equal(forecast(chained, h = 5)$.mean, forecast(one_shot, h = 5)$.mean)
})

test_that("stream.fbl_ARFIMA validates that new_data starts immediately after the trained data", {
  skip_if_not_installed("fracdiff")
  fit <- head(sunspot_tbl, -20) %>% model(arfima = ARFIMA(value ~ pdq(1, NULL, 0)))
  expect_error(stream(fit, tail(sunspot_tbl, 19)), "must start one step beyond")
})
