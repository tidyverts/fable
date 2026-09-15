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

test_that("ETS with missing values", {
  UK_missing <- UKLungDeaths
  UK_missing[["mdeaths"]][3:5] <- NA
  fit_missing <- expect_no_error(UK_missing |> model(ETS(mdeaths)))
  # A seasonal model should still be selected despite the missing values
  expect_true(fit_missing[[1]][[1]]$fit$spec$seasontype != "N")

  USAccDeaths_miss <- USAccDeaths_tbl
  USAccDeaths_miss$value[c(10, 14, 15)] <- NA
  fable_fit <- USAccDeaths_miss |> model(ets = ETS(value))

  USAccDeaths_miss <- fable_fit |>
    interpolate(USAccDeaths_miss)
  expect_false(
    any(is.na(USAccDeaths_miss$value))
  )
  expect_equal(
    USAccDeaths_tbl$value[-c(10, 14, 15)],
    USAccDeaths_miss$value[-c(10, 14, 15)]
  )
})


# stream() should exactly match refitting the model (with the same smoothing
# parameters and initial states) to the complete series.
expect_stream_matches_refit <- function(streamed, refitted) {
  s <- streamed[[1]][[1]]$fit
  r <- refitted[[1]][[1]]$fit
  expect_equal(tidy(streamed)$estimate, tidy(refitted)$estimate)
  expect_equal(fitted(streamed)[[".fitted"]], fitted(refitted)[[".fitted"]])
  expect_equal(residuals(streamed)[[".resid"]], residuals(refitted)[[".resid"]])
  expect_equal(s$states, r$states)
  expect_equal(glance(streamed), glance(refitted))
  expect_equal(s$amse, r$amse)
}

test_that("stream.ETS matches a full refit", {
  tr <- USAccDeaths_tbl %>% head(-12)
  nw <- USAccDeaths_tbl %>% tail(12)

  for (spec in list(
    value ~ error("A") + trend("N") + season("A"),
    value ~ error("A") + trend("Ad") + season("A"),
    value ~ error("M") + trend("A") + season("M"),
    value ~ error("M") + trend("Ad") + season("M"),
    value ~ error("A") + trend("A") + season("N"),
    value ~ error("M") + trend("N") + season("N")
  )) {
    fit <- tr %>% model(ets = ETS(!!spec))
    streamed <- fit %>% stream(nw)
    refitted <- fit %>% refit(USAccDeaths_tbl, reinitialise = FALSE)

    expect_stream_matches_refit(streamed, refitted)
    expect_s3_class(streamed$ets[[1]]$fit$states, "tbl_ts")
    expect_equal(NROW(augment(streamed)), NROW(USAccDeaths_tbl))
  }
})

test_that("stream.ETS chained calls match a single stream()", {
  fit <- USAccDeaths_tbl %>%
    head(-12) %>%
    model(ets = ETS(value ~ error("A") + trend("Ad") + season("A")))
  refitted <- fit %>% refit(USAccDeaths_tbl, reinitialise = FALSE)

  # Streaming one observation at a time exercises the multi-step (AMSE)
  # lookback across every boundary.
  streamed <- fit
  for (i in seq_len(12)) {
    streamed <- streamed %>% stream(USAccDeaths_tbl[NROW(USAccDeaths_tbl) - 12 + i, ])
  }
  expect_stream_matches_refit(streamed, refitted)

  # Forecasts continue from the streamed states
  expect_equal(
    forecast(streamed, h = 6)$.mean,
    forecast(refitted, h = 6)$.mean
  )
})

test_that("stream.ETS handles missing values", {
  dat <- USAccDeaths_tbl
  dat$value[c(20, 65, 70)] <- NA
  fit <- dat %>%
    head(-12) %>%
    model(ets = ETS(value ~ error("A") + trend("A") + season("A")))
  streamed <- fit %>% stream(dat %>% tail(12) %>% head(6)) %>% stream(dat %>% tail(6))
  refitted <- fit %>% refit(dat, reinitialise = FALSE)

  expect_stream_matches_refit(streamed, refitted)
})

test_that("stream.ETS must start immediately after the trained data", {
  fit <- USAccDeaths_tbl %>%
    head(-12) %>%
    model(ets = ETS(value ~ error("A") + trend("N") + season("A")))
  expect_error(
    fit %>% stream(USAccDeaths_tbl %>% tail(11)),
    "must start one step beyond the end of"
  )
})
