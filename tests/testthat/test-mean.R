context("test-mean")

test_that("MEAN", {
  fit <- USAccDeaths_tbl %>% model(mean = MEAN(value))
  expect_identical(
    fitted(fit) %>% select(index, .fitted),
    USAccDeaths_tbl %>% transmute(.fitted = mean(value))
  )

  expect_equivalent(
    glance(fit)$sigma2, var(scale(USAccDeaths_tbl$value, scale = FALSE))
  )

  expect_identical(
    residuals(fit)[[".resid"]], USAccDeaths_tbl$value - mean(USAccDeaths_tbl$value)
  )

  expect_identical(
    tidy(fit)$estimate, mean(USAccDeaths_tbl$value)
  )

  expect_output(report(fit), "Mean: 8788\\.7917")

  fc <- fit %>% forecast(h = 3)

  expect_identical(
    fc_mean(fc$value),
    rep(mean(USAccDeaths_tbl$value), 3)
  )

  fc_sim <- fit %>% forecast(h = 3, bootstrap = TRUE, times = 5)

  # expect_identical(
  #   fc$value, fc_sim$value
  # )
})

test_that("stream.model_mean", {
  library(tsibble)
  lung_deaths_male <- as_tsibble(mdeaths)

  train <- lung_deaths_male %>% filter(index < yearmonth("1979 Jan"))
  new_obs <- lung_deaths_male %>% filter(index >= yearmonth("1979 Jan"))

  # Fixed mean: streaming should exactly match a fixed-mean refit
  fit <- train %>% model(mean = MEAN(value))
  streamed <- fit %>% stream(new_obs)
  refit_fit <- fit %>% refit(lung_deaths_male, reestimate = FALSE)

  expect_equal(tidy(streamed)$estimate, tidy(fit)$estimate)
  expect_equal(
    fitted(streamed)[[".fitted"]],
    fitted(refit_fit)[[".fitted"]]
  )
  expect_equal(
    residuals(streamed)[[".resid"]],
    residuals(refit_fit)[[".resid"]]
  )
  expect_equal(glance(streamed)$sigma2, glance(refit_fit)$sigma2)
  expect_equal(tidy(streamed), tidy(refit_fit))
  expect_equal(
    streamed %>% forecast(h = 6),
    refit_fit %>% forecast(h = 6)
  )

  # Rolling window mean: streaming should exactly match estimating on the full series
  fit_w <- train %>% model(mean = MEAN(value ~ window(size = 6)))
  streamed_w <- fit_w %>% stream(new_obs)
  refit_w <- fit_w %>% refit(lung_deaths_male, reestimate = TRUE)

  expect_equal(
    fitted(streamed_w)[[".fitted"]],
    fitted(refit_w)[[".fitted"]]
  )
  expect_equal(
    residuals(streamed_w)[[".resid"]],
    residuals(refit_w)[[".resid"]]
  )
  expect_equal(glance(streamed_w)$sigma2, glance(refit_w)$sigma2)
  expect_equal(tidy(streamed_w), tidy(refit_w))
  expect_equal(
    streamed_w %>% forecast(h = 6),
    refit_w %>% forecast(h = 6)
  )

  # Rolling window larger than the trained data should still match
  fit_short <- lung_deaths_male %>%
    filter(index < yearmonth("1974 Apr")) %>%
    model(mean = MEAN(value ~ window(size = 6)))
  streamed_short <- fit_short %>%
    stream(lung_deaths_male %>% filter(index >= yearmonth("1974 Apr")))
  expect_equal(
    fitted(streamed_short)[[".fitted"]],
    fitted(refit_w)[[".fitted"]]
  )
  expect_equal(tidy(streamed_short), tidy(refit_w))

  # Chained stream() calls should give the same result as one large stream()
  mid <- lung_deaths_male %>% filter(index < yearmonth("1978 Jul"))
  part1 <- lung_deaths_male %>%
    filter(index >= yearmonth("1978 Jul"), index < yearmonth("1979 Jan"))
  part2 <- new_obs

  chained <- mid %>%
    model(mean = MEAN(value)) %>%
    stream(part1) %>%
    stream(part2)
  full <- mid %>%
    model(mean = MEAN(value)) %>%
    stream(dplyr::bind_rows(part1, part2))
  expect_equal(fitted(chained), fitted(full))
  expect_equal(residuals(chained), residuals(full))
  expect_equal(glance(chained), glance(full))

  chained_w <- fit_w
  for (i in seq_len(NROW(new_obs))) {
    chained_w <- chained_w %>% stream(new_obs[i, ])
  }
  expect_equal(fitted(chained_w), fitted(streamed_w))
  expect_equal(residuals(chained_w), residuals(streamed_w))
  expect_equal(glance(chained_w), glance(streamed_w))
  expect_equal(
    chained_w %>% forecast(h = 6),
    streamed_w %>% forecast(h = 6)
  )

  # Streaming must start immediately after the trained data
  expect_error(
    fit %>% stream(lung_deaths_male %>% filter(index >= yearmonth("1979 Feb"))),
    "must start one step beyond the end of"
  )
})
