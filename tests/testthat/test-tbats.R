context("test-tbats.R")

test_that("Automatic TBATS selection", {
  fit <- USAccDeaths_tbl %>% model(tbats = TBATS(value))

  expect_s3_class(fit$tbats[[1]]$fit, "TBATS")
  expect_identical(fit$tbats[[1]]$fit$spec$period[[1]], 12)

  # Coefficients are sensible (admissible smoothing parameter)
  par <- tidy(fit)
  expect_true(all(c("alpha", "gamma1_1", "gamma2_1") %in% par$term))

  fc <- fit %>% forecast(h = 24)
  expect_equal(NROW(fc), 24)
  # Prediction intervals should (weakly) widen with the forecast horizon
  sd <- distributional::variance(fc$value) %>% sqrt()
  expect_true(all(diff(sd) >= -1e-8))
})

test_that("Manual TBATS specification", {
  fit <- USAccDeaths_tbl %>%
    model(tbats = TBATS(value ~ trend("Ad") + season(period = 12, harmonics = 3) + pdq(p = 0, q = 0)))

  spec <- fit$tbats[[1]]$fit$spec
  expect_true(spec$use_trend)
  expect_true(spec$damped)
  expect_identical(spec$k[[1]], 3L)
  expect_identical(spec$p, 0L)
  expect_identical(spec$q, 0L)

  expect_identical(
    model_sum(fit$tbats[[1]]),
    sprintf(
      "TBATS({0,0}, %s, {<12,3>})",
      format(round(tidy(fit)$estimate[tidy(fit)$term == "phi"], 3))
    )
  )
})

test_that("TBATS disables ARMA errors with pdq(0, 0)", {
  fit <- USAccDeaths_tbl %>% model(tbats = TBATS(value ~ pdq(p = 0, q = 0)))
  expect_identical(fit$tbats[[1]]$fit$spec$p, 0L)
  expect_identical(fit$tbats[[1]]$fit$spec$q, 0L)
})

test_that("TBATS generate() is reproducible and refit() behaves correctly", {
  fit <- USAccDeaths_tbl %>% model(tbats = TBATS(value))

  set.seed(1)
  sim1 <- generate(fit, h = 12, times = 1)
  set.seed(1)
  sim2 <- generate(fit, h = 12, times = 1)
  expect_equal(sim1$.sim, sim2$.sim)

  # reinitialise = TRUE re-estimates parameters (won't generally match)
  refit_new <- refit(fit, USAccDeaths_tbl, reinitialise = TRUE)
  expect_s3_class(refit_new$tbats[[1]]$fit, "TBATS")

  # reinitialise = FALSE keeps the same structure/coefficients when applied
  # to the same training data - fitted values should match closely.
  refit_same <- refit(fit, USAccDeaths_tbl, reinitialise = FALSE)
  expect_identical(tidy(fit)$estimate, tidy(refit_same)$estimate)
  expect_equal(fitted(fit)$.fitted, fitted(refit_same)$.fitted, tolerance = 1e-6)
})

test_that("TBATS components reconstruct the response", {
  fit <- USAccDeaths_tbl %>% model(tbats = TBATS(value))
  cmp <- components(fit)

  # State rows are the *post-update* state, so reconstructing y[t] needs the
  # previous row's level/season (see components.TBATS()'s `aliases` formula).
  n <- NROW(cmp)
  reconstructed <- cmp$level[-n] + cmp$season1[-n] + cmp$remainder[-1]
  actual <- cmp$value[-1]
  expect_equal(reconstructed, actual, tolerance = 1e-6)
})

test_that("TBATS works with a box_cox() transformed response", {
  fit <- USAccDeaths_tbl %>% model(tbats = TBATS(box_cox(value, 0.3)))
  expect_s3_class(fit$tbats[[1]]$fit, "TBATS")

  fc <- fit %>% forecast(h = 3)
  # Back-transformed forecasts should be on the original (untransformed) scale
  expect_true(all(fc_mean(fc$value) > 1000))
})

test_that("stream.TBATS matches a full refit closely", {
  full_data <- UKLungDeaths
  early <- full_data %>% dplyr::filter(index < yearmonth("1979 Jan"))
  late <- full_data %>% dplyr::filter(index >= yearmonth("1979 Jan"))

  fit <- early %>% model(tbats = TBATS(mdeaths))
  streamed <- fit %>% stream(late)

  expect_equal(NROW(streamed$tbats[[1]]$fit$est), NROW(full_data))
  # Coefficients are unchanged by streaming (no re-estimation)
  expect_identical(tidy(fit)$estimate, tidy(streamed)$estimate)

  full_fit <- full_data %>% model(tbats = TBATS(mdeaths))
  expect_equal(glance(streamed)$AIC, glance(full_fit)$AIC, tolerance = 0.05)
})

test_that("TBATS forecasts do not depend on the order of seasonal periods", {
  fc_asc <- USAccDeaths_tbl %>%
    model(TBATS(value ~ season(period = c(4, 12)) + pdq(0, 0))) %>%
    forecast(h = 6)
  fc_desc <- USAccDeaths_tbl %>%
    model(TBATS(value ~ season(period = c(12, 4)) + pdq(0, 0))) %>%
    forecast(h = 6)
  expect_equal(fc_asc$.mean, fc_desc$.mean)
  expect_equal(distributional::variance(fc_asc$value), distributional::variance(fc_desc$value))
})
