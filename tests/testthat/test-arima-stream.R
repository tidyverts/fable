context("test-arima-stream.R")

test_that("stream.ARIMA matches a full refit (no differencing)", {
  tr <- USAccDeaths_tbl %>% head(-12)
  nw <- USAccDeaths_tbl %>% tail(12)

  fit <- tr %>% model(arima = ARIMA(value ~ 1 + pdq(2, 0, 1) + PDQ(0, 0, 0)))
  streamed <- fit %>% stream(nw)

  coef_fixed <- fit$arima[[1]]$fit$model$coef
  intercept <- fable:::arima_constant(NROW(USAccDeaths_tbl), 0, 0, 1)
  stats_fit <- arima(USAccDeaths, order = c(2, 0, 1),
    xreg = matrix(intercept, dimnames = list(NULL, "constant")),
    include.mean = FALSE, fixed = coef_fixed, method = "ML"
  )

  # Coefficients are unchanged by streaming (no re-estimation)
  expect_identical(tidy(fit)$estimate, tidy(streamed)$estimate)

  # Fitted values, residuals and fit statistics exactly match a full refit
  # with the same (fixed) coefficients, since there's no differencing.
  expect_equivalent(fitted(streamed)$.fitted, unclass(USAccDeaths - residuals(stats_fit)))
  expect_equivalent(residuals(streamed)$.resid, unclass(residuals(stats_fit)))
  expect_equal(glance(streamed)$log_lik, stats_fit$loglik)

  npar <- length(coef_fixed) + 1
  nstar <- length(USAccDeaths)
  sigma2 <- sum(residuals(stats_fit)^2) / (nstar - npar + 1)
  expect_equal(glance(streamed)$sigma2, sigma2)
  expect_equal(glance(streamed)$AIC, -2 * stats_fit$loglik + 2 * npar)

  # object$model (the underlying Arima object) is kept in sync too, since
  # predict.Arima() relies on `sigma2` for forecast standard errors.
  expect_equal(streamed$arima[[1]]$fit$model$sigma2, sigma2)
  expect_equal(streamed$arima[[1]]$fit$model$loglik, stats_fit$loglik)
})

test_that("stream.ARIMA fitted values match a full refit (with differencing)", {
  tr <- USAccDeaths_tbl %>% head(-12)
  nw <- USAccDeaths_tbl %>% tail(12)

  fit <- tr %>% model(arima = ARIMA(value ~ pdq(0, 1, 1) + PDQ(0, 1, 1)))
  streamed <- fit %>% stream(nw)

  coef_fixed <- fit$arima[[1]]$fit$model$coef
  stats_fit <- arima(USAccDeaths, order = c(0, 1, 1), seasonal = list(order = c(0, 1, 1), period = 12),
    fixed = coef_fixed, method = "ML"
  )

  # .fitted/.resid/sigma2 remain exact even with (seasonal) differencing, as
  # these come from the incrementally-updated Kalman filter state rather
  # than the (approximate, see ?stream.ARIMA) combined log-likelihood.
  expect_equivalent(fitted(streamed)$.fitted, unclass(USAccDeaths - residuals(stats_fit)))
  expect_equivalent(residuals(streamed)$.resid, unclass(residuals(stats_fit)))

  npar <- length(coef_fixed) + 1
  nstar <- length(USAccDeaths) - 1 - 12
  sigma2 <- sum(residuals(stats_fit)^2, na.rm = TRUE) / (nstar - npar + 1)
  expect_equal(glance(streamed)$sigma2, sigma2)

  # log_lik/AIC/AICc/BIC are only approximate once differencing is used, but
  # should still be reasonably close to the exact value.
  expect_equal(glance(streamed)$log_lik, stats_fit$loglik, tolerance = 0.05)
})

test_that("stream.ARIMA handles exogenous regressors and rank-deficient xreg", {
  tr <- UKLungDeaths %>% head(-12)
  nw <- UKLungDeaths %>% tail(12)

  fit <- tr %>% model(model = ARIMA(mdeaths ~ 1 + fdeaths + PDQ(P = 0, Q = 0)))
  streamed <- tryCatch(stream(fit, nw), error = function(e) e)
  expect_false(inherits(streamed, "error"))
  expect_identical(tidy(fit)$estimate, tidy(streamed)$estimate)
  expect_equal(nrow(streamed$model[[1]]$fit$est), nrow(UKLungDeaths))

  # A rank-deficient regressor dropped at training time should not break
  # streaming (the same columns must be dropped consistently).
  UKLungDeaths_rd <- UKLungDeaths
  UKLungDeaths_rd$fdeaths2 <- UKLungDeaths_rd$fdeaths * 2
  tr_rd <- UKLungDeaths_rd %>% head(-12)
  nw_rd <- UKLungDeaths_rd %>% tail(12)
  fit_rd <- suppressWarnings(
    tr_rd %>% model(model = ARIMA(mdeaths ~ 1 + fdeaths + fdeaths2 + pdq(1, 0, 0) + PDQ(P = 0, Q = 0)))
  )
  streamed_rd <- tryCatch(stream(fit_rd, nw_rd), error = function(e) e)
  expect_false(inherits(streamed_rd, "error"))
})

test_that("stream.ARIMA validates that new_data starts immediately after the trained data", {
  tr <- USAccDeaths_tbl %>% head(-12)
  fit <- tr %>% model(arima = ARIMA(value ~ 1 + pdq(1, 0, 0) + PDQ(0, 0, 0)))

  gap <- USAccDeaths_tbl %>% tail(11) # skips one observation
  expect_error(
    stream(fit, gap),
    "must start one step beyond"
  )

  contiguous <- USAccDeaths_tbl %>% tail(12)
  expect_error(stream(fit, contiguous), NA)
})

test_that("Chained streaming matches one-shot streaming", {
  tr <- USAccDeaths_tbl %>% head(-24)
  fit <- tr %>% model(arima = ARIMA(value ~ 1 + pdq(1, 0, 1) + PDQ(0, 0, 0)))

  chained <- fit %>%
    stream(USAccDeaths_tbl %>% tail(24) %>% head(12)) %>%
    stream(USAccDeaths_tbl %>% tail(12))
  one_shot <- fit %>% stream(USAccDeaths_tbl %>% tail(24))

  expect_equal(fitted(chained)$.fitted, fitted(one_shot)$.fitted)
  expect_equal(residuals(chained)$.resid, residuals(one_shot)$.resid)
  expect_equal(glance(chained)$log_lik, glance(one_shot)$log_lik)
  expect_equal(glance(chained)$AIC, glance(one_shot)$AIC)
})

test_that("stream.ARIMA updates forecast standard errors via sigma2", {
  tr <- USAccDeaths_tbl %>% head(-24)
  fit <- tr %>% model(arima = ARIMA(value ~ 1 + pdq(1, 0, 0) + PDQ(0, 0, 0)))
  pre_sigma2 <- fit$arima[[1]]$fit$model$sigma2

  streamed <- fit %>% stream(USAccDeaths_tbl %>% tail(24))
  post_sigma2 <- streamed$arima[[1]]$fit$model$sigma2

  expect_false(isTRUE(all.equal(pre_sigma2, post_sigma2)))

  fc <- forecast(streamed, h = 2)
  expect_equal(
    sqrt(distributional::variance(fc$value)),
    unname(sqrt(stats::KalmanForecast(2, streamed$arima[[1]]$fit$model$model)$var * post_sigma2))
  )
})

test_that("stream.ARIMA resolves lag() regressors on a single stream() call", {
  tr <- UKLungDeaths %>% head(-24)
  fit <- tr %>% model(model = ARIMA(mdeaths ~ 1 + lag(fdeaths) + pdq(1, 0, 0) + PDQ(P = 0, Q = 0)))

  # `refit(reestimate = FALSE)` reapplies the trained coefficients to the
  # complete series in one pass, giving the correct ground truth to compare
  # against (it always sees the full history, so isn't affected by the
  # streaming boundary that `lag()` needs `self$recent_data` for).
  ground_truth <- refit(fit, UKLungDeaths, reestimate = FALSE)

  one_shot <- fit %>% stream(UKLungDeaths %>% tail(24))
  expect_equal(residuals(one_shot)$.resid, residuals(ground_truth)$.resid)
  # Only the very first observation of the whole series is unresolvable -
  # there's no earlier data for lag() to use, streamed or otherwise.
  expect_identical(which(is.na(residuals(one_shot)$.resid)), 1L)

  # Chained (repeated) streaming resolves lag() correctly too.
  chained <- fit %>%
    stream(UKLungDeaths %>% tail(24) %>% head(12)) %>%
    stream(UKLungDeaths %>% tail(12))
  expect_equal(residuals(chained)$.resid, residuals(ground_truth)$.resid)
})
