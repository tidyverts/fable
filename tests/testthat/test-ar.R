context("test-ar.R")

lh_tbl <- as_tsibble(lh)
lh_train <- lh_tbl %>% filter(index <= 40)
lh_new <- lh_tbl %>% filter(index > 40)

# Streaming new_obs onto fit should match a fixed-coefficient refit to full_data
expect_stream_matches_refit <- function(fit, new_obs, full_data) {
  streamed <- fit %>% stream(new_obs)
  refitted <- fit %>% refit(full_data, reestimate = FALSE)

  expect_equal(tidy(streamed)$estimate, tidy(fit)$estimate)
  expect_equal(tidy(streamed)$std.error, tidy(fit)$std.error)
  expect_equal(tidy(streamed), tidy(refitted))
  expect_equal(fitted(streamed)[[".fitted"]], fitted(refitted)[[".fitted"]])
  expect_equal(residuals(streamed)[[".resid"]], residuals(refitted)[[".resid"]])
  expect_equal(
    residuals(streamed, type = "regression")[[".resid"]],
    residuals(refitted, type = "regression")[[".resid"]]
  )
  expect_equal(glance(streamed)$sigma2, glance(refitted)$sigma2)
  expect_equal(glance(streamed)$dof, glance(refitted)$dof)

  # refit() treats all coefficients as fixed, so its IC has no parameter penalty
  npar <- fit[[1]][[1]]$fit$npar
  n <- nrow(full_data)
  expect_equal(glance(streamed)$AIC, glance(refitted)$AIC + 2 * npar)
  expect_equal(glance(streamed)$BIC, glance(refitted)$BIC + npar * log(n))
  expect_equal(
    glance(streamed)$AICc,
    glance(refitted)$AICc + 2 * npar + 2 * npar * (npar + 1) / (n - npar - 1)
  )

  fc_streamed <- forecast(streamed, h = 6)
  fc_refitted <- forecast(refitted, h = 6)
  resp <- response_vars(fc_streamed)
  expect_equal(fc_mean(fc_streamed[[resp]]), fc_mean(fc_refitted[[resp]]))
  expect_equal(
    distributional::variance(fc_streamed[[resp]]),
    distributional::variance(fc_refitted[[resp]])
  )

  invisible(streamed)
}

test_that("stream.AR with constant", {
  fit <- lh_train %>% model(ar = AR(value ~ order(3)))
  expect_equal(fit$ar[[1]]$fit$npar, 4)
  streamed <- expect_stream_matches_refit(fit, lh_new, lh_tbl)
  expect_equal(nrow(fitted(streamed)), nrow(lh_tbl))
})

test_that("stream.AR without constant", {
  fit <- lh_train %>% model(ar = AR(value ~ 0 + order(2)))
  expect_equal(fit$ar[[1]]$fit$npar, 2)
  expect_stream_matches_refit(fit, lh_new, lh_tbl)
})

test_that("stream.AR with exogenous regressors", {
  fit <- UKLungDeaths %>%
    head(-12) %>%
    model(ar = AR(mdeaths ~ fdeaths + order(2)))
  expect_stream_matches_refit(fit, tail(UKLungDeaths, 12), UKLungDeaths)

  fit_trend <- lh_train %>% model(ar = AR(value ~ trend() + order(2)))
  expect_stream_matches_refit(fit_trend, lh_new, lh_tbl)
})

test_that("Chained stream.AR calls match one large stream", {
  fit <- lh_tbl %>%
    filter(index <= 36) %>%
    model(ar = AR(value ~ order(3)))

  chained <- fit %>%
    stream(lh_tbl %>% filter(index > 36, index <= 40)) %>%
    stream(lh_new)
  one_shot <- fit %>% stream(lh_tbl %>% filter(index > 36))

  expect_equal(fitted(chained), fitted(one_shot))
  expect_equal(residuals(chained), residuals(one_shot))
  expect_equal(glance(chained), glance(one_shot))
  expect_equal(fit %>% refit(lh_tbl, reestimate = FALSE) %>% glance() %>% .$sigma2, glance(chained)$sigma2)

  # Single observation streams
  single <- fit
  for (i in 37:48) single <- single %>% stream(lh_tbl %>% filter(index == i))
  expect_equal(fitted(single), fitted(one_shot))
  expect_equal(glance(single), glance(one_shot))
  expect_equal(
    fc_mean(forecast(single, h = 3)$value),
    fc_mean(forecast(one_shot, h = 3)$value)
  )
})

test_that("stream.AR validates that new_data starts immediately after the trained data", {
  fit <- lh_train %>% model(ar = AR(value ~ order(3)))
  expect_error(
    fit %>% stream(lh_tbl %>% filter(index > 41)),
    "must start one step beyond the end of"
  )
  expect_error(fit %>% stream(lh_new), NA)
})

test_that("Fixed AR regression coefficients are in the original scale", {
  fit <- lh_train %>% model(ar = AR(value ~ order(3)))
  fixed <- lh_train %>%
    model(ar = AR(value ~ order(3, fixed = list(ar1 = 0.5)) + xreg(fixed = list(constant = 1))))
  expect_equal(tidy(fixed)$estimate[1:2], c(1, 0.5))
  expect_equal(
    tidy(fit %>% refit(lh_tbl, reestimate = FALSE))$estimate,
    tidy(fit)$estimate
  )
})
