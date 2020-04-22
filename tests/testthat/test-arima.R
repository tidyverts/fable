context("test-arima.R")

stats_fit <- arima(USAccDeaths, c(0, 1, 1), list(order = c(0, 1, 1), 12))
test_that("Automatic ARIMA selection", {
  skip_if_not_installed("feasts")
  # Automatic model selection
  fable_fit <- USAccDeaths_tbl %>%
    model(arima = ARIMA(value ~ pdq(0:1, 0:1, 0:1) + PDQ(0:1, 0:1, 0:1)))

  expect_identical(
    coef(fable_fit$arima[[1]]$fit$model),
    coef(stats_fit)
  )

  # Automatic (approximate) model selection
  fable_fit_approx <- USAccDeaths_tbl %>%
    model(arima = ARIMA(value ~ pdq(0:1, 0:1, 0:1) + PDQ(0:1, 0:1, 0:1), approximation = TRUE))

  expect_output(
    report(fable_fit_approx),
    "ARIMA\\(0,1,1\\)\\(1,1,0\\)\\[12\\]"
  )

  # Automatic d/D selection
  fit <- USAccDeaths_tbl %>% model(ARIMA(value ~ pdq(p = 0, q = 1) + PDQ(P = 0, Q = 1)))
  expect_identical(
    model_sum(fit[[1]][[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )
})

test_that("Manual ARIMA selection", {
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(model = ARIMA(value ~ 0 + pdq(0, 1, 1) + PDQ(0, 1, 1)))

  expect_identical(
    coef(fable_fit$model[[1]]$fit$model),
    coef(stats_fit)
  )

  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "ARIMA(0,1,1)(0,1,1)[12]"
  )

  fable_fc <- fable_fit %>% forecast()
  stats_fc <- stats_fit %>% predict(24)

  expect_equivalent(
    fc_mean(fable_fc$value),
    unclass(stats_fc$pred)
  )

  expect_equivalent(
    tidy(fable_fit)$estimate, coef(stats_fit)
  )

  expect_equivalent(
    fitted(fable_fit)$.fitted, unclass(USAccDeaths - residuals(stats_fit))
  )

  expect_equivalent(
    residuals(fable_fit)$.resid, unclass(residuals(stats_fit))
  )

  expect_equivalent(
    fitted(fable_fit),
    fitted(refit(fable_fit, USAccDeaths_tbl))
  )

  expect_output(
    report(fable_fit),
    "log likelihood=-425.44"
  )

  USAccDeaths_miss <- USAccDeaths_tbl
  USAccDeaths_miss$value[c(10, 14, 15)] <- NA
  USAccDeaths_miss <- fable_fit %>%
    interpolate(USAccDeaths_miss)
  expect_false(
    any(is.na(USAccDeaths_miss$value))
  )
  expect_equal(
    USAccDeaths_tbl$value[-c(10, 14, 15)],
    USAccDeaths_miss$value[-c(10, 14, 15)]
  )
})

test_that("Fixed ARIMA coefficients", {
  # Manual model selection
  fable_fit <- USAccDeaths_tbl %>% model(model = ARIMA(value ~ xreg(1, fixed = list(constant = 20)) + pdq(0, 0, 1, fixed = list(ma1 = 0.3)) + PDQ(0, 1, 1, fixed = list(sma1 = 3))))
  
  expect_identical(
    tidy(fable_fit)$estimate,
    c(0.3, 3, 20)
  )
})

test_that("ARIMA with bad inputs", {
  expect_warning(
    UKLungDeaths %>%
      model(ARIMA(mdeaths ~ 1 + pdq(2, 0, 0) + pdq(d = 1) + PDQ(2, 1, 0))),
    "Only one special for `pdq\\(\\)` and `PDQ\\(\\)` is allowed"
  )

  expect_warning(
    UKLungDeaths %>%
      model(ARIMA(vars(mdeaths, fdeaths))),
    "Only univariate responses are supported by ARIMA"
  )

  expect_warning(
    UKLungDeaths %>%
      model(ARIMA(mdeaths ~ pdq(2, 2, 0) + PDQ(2, 1, 0))),
    "Having 3 or more differencing operations is not recommended"
  )

  expect_warning(
    UKLungDeaths %>%
      model(ARIMA(mdeaths ~ pdq(2, 0, 0) + PDQ(2, 2, 0))),
    "Having more than one seasonal differences is not recommended"
  )

  expect_warning(
    fit <- UKLungDeaths %>%
      model(poly = ARIMA(mdeaths ~ 1 + pdq(2, 1, 0) + PDQ(2, 1, 0))),
    "Model specification induces a quadratic or higher order polynomial trend"
  )

  expect_identical(
    model_sum(fit$poly[[1]]),
    "ARIMA(2,1,0)(2,1,0)[12] w/ poly"
  )
})

test_that("ARIMA with xregs", {
  skip_if_not_installed("feasts")
  tr <- UKLungDeaths %>% head(-12)
  ts <- UKLungDeaths %>% tail(12)
  fable_fit <- tr %>% model(model = ARIMA(mdeaths ~ 1 + fdeaths + PDQ(P = 0, Q = 0)))
  stats_fit <- arima(head(mdeaths, -12), c(1, 1, 1),
    xreg = data.frame(fdeaths = head(fdeaths, -12), intercept = seq_len(60))
  )

  expect_equivalent(
    coef(fable_fit$model[[1]]$fit$model),
    coef(stats_fit)
  )

  fable_fc <- fable_fit %>% forecast(ts)
  stats_fc <- stats_fit %>% predict(12,
    newxreg = data.frame(fdeaths = tail(fdeaths, 12), intercept = 61:72)
  )

  expect_equivalent(
    fc_mean(fable_fc$mdeaths),
    unclass(stats_fc$pred)
  )

  expect_identical(
    model_sum(fable_fit$model[[1]]),
    "LM w/ ARIMA(1,1,1) errors"
  )

  fable_fit <- tr %>%
    model(model = ARIMA(mdeaths ~ 1 + lag(fdeaths) + PDQ(P = 0, Q = 0)))
  expect_equal(
    model_sum(fable_fit$model[[1]]),
    "LM w/ ARIMA(2,0,1) errors"
  )
  fable_fc <- fable_fit %>%
    forecast(ts)
  expect_true(
    !any(is.na(fable_fc$mdeaths))
  )
})
