context("test-theta.R")

# Manually apply a Theta model's fixed parameters to new observations
theta_manual_stream <- function(mdl, y) {
  h <- length(y)
  m <- mdl$period
  additive <- mdl$dcmp == "additive"
  seas <- if (m > 1L) rep(mdl$season, length.out = h) else rep(1, h)
  y_sa <- if (m == 1L) y else if (additive) y - seas else y / seas

  l <- unname(mdl$lT)
  fitted_sa <- numeric(h)
  for (i in seq_len(h)) {
    fitted_sa[i] <- l
    if (!is.na(y_sa[i])) l <- l + unname(mdl$alpha) * (y_sa[i] - l)
  }
  fitted <- if (m == 1L) fitted_sa else if (additive) fitted_sa + seas else fitted_sa * seas

  n_old <- length(mdl$resid)
  list(
    fitted = fitted,
    resid = y - fitted,
    lT = l,
    sigma2 = (mdl$sigma2 * (n_old - 2) + sum((y_sa - fitted_sa)^2, na.rm = TRUE)) / (n_old + h - 2),
    season = if (m > 1L) rep(mdl$season, length.out = h + m)[h + seq_len(m)] else NULL
  )
}

test_that("stream.fable_theta matches a manual SES recursion (non-seasonal)", {
  lh_tbl <- as_tsibble(lh)
  tr <- head(lh_tbl, -8)
  nw <- tail(lh_tbl, 8)

  fit <- tr %>% model(theta = THETA(value))
  mdl <- fit$theta[[1]]$fit
  expect_equal(mdl$period, 1L)

  streamed <- fit %>% stream(nw)
  smdl <- streamed$theta[[1]]$fit
  expected <- theta_manual_stream(mdl, nw$value)

  # Parameters are unchanged
  expect_equal(tidy(streamed)$estimate, tidy(fit)$estimate)

  n_old <- NROW(tr)
  expect_equal(fitted(streamed)[[".fitted"]][-seq_len(n_old)], expected$fitted)
  expect_equal(residuals(streamed)[[".resid"]][-seq_len(n_old)], expected$resid)
  expect_equal(fitted(streamed)[[".fitted"]][seq_len(n_old)], fitted(fit)[[".fitted"]])
  expect_equal(smdl$lT, expected$lT)
  expect_equal(NROW(augment(streamed)), NROW(lh_tbl))

  # sigma2 uses the same formula as training over the complete residuals
  expect_equal(glance(streamed)$sigma2, expected$sigma2)
  expect_equal(
    glance(streamed)$sigma2,
    sum(residuals(streamed)[[".resid"]]^2) / (NROW(lh_tbl) - 2)
  )

  # Forecasts continue from the streamed level and use the full series length
  fc <- forecast(streamed, h = 4)
  alpha <- unname(mdl$alpha)
  n <- NROW(lh_tbl)
  expect_equal(
    fc$.mean,
    expected$lT + unname(mdl$drift) * (0:3 + (1 - (1 - alpha)^n) / alpha)
  )
  expect_equal(NROW(fc), 4)
  expect_equal(fc$index, max(lh_tbl$index) + 1:4)
})

test_that("stream.fable_theta matches a manual SES recursion (seasonal)", {
  tr <- USAccDeaths_tbl %>% head(-12)
  nw <- USAccDeaths_tbl %>% tail(12)

  for (method in c("multiplicative", "additive")) {
    fit <- tr %>% model(theta = THETA(value ~ season(method = !!method)))
    mdl <- fit$theta[[1]]$fit
    expect_equal(unname(mdl$period), 12L)
    expect_equal(mdl$dcmp, method)

    streamed <- fit %>% stream(nw)
    smdl <- streamed$theta[[1]]$fit
    expected <- theta_manual_stream(mdl, nw$value)

    expect_equal(tidy(streamed)$estimate, tidy(fit)$estimate)
    expect_equal(smdl$season, mdl$season)

    n_old <- NROW(tr)
    expect_equal(fitted(streamed)[[".fitted"]][-seq_len(n_old)], expected$fitted)
    expect_equal(residuals(streamed)[[".resid"]][-seq_len(n_old)], expected$resid)
    expect_equal(smdl$lT, expected$lT)
    expect_equal(glance(streamed)$sigma2, expected$sigma2)
    expect_equal(NROW(augment(streamed)), NROW(USAccDeaths_tbl))

    if (method == "additive") {
      # Additive residuals are unaffected by reseasonalisation
      expect_equal(
        glance(streamed)$sigma2,
        sum(residuals(streamed)[[".resid"]]^2) / (NROW(USAccDeaths_tbl) - 2)
      )
    }

    # Forecasts continue from the streamed level and reseasonalise
    h <- 18
    fc <- forecast(streamed, h = h)
    alpha <- unname(mdl$alpha)
    n <- NROW(USAccDeaths_tbl)
    fc_sa <- expected$lT + unname(mdl$drift) * (0:(h - 1) + (1 - (1 - alpha)^n) / alpha)
    seas <- rep(expected$season, length.out = h)
    expect_equal(fc$.mean, if (method == "additive") fc_sa + seas else fc_sa * seas)
    expect_equal(NROW(fc), h)
    expect_equal(fc$index, max(USAccDeaths_tbl$index) + 1:h)
  }
})

test_that("stream.fable_theta chained calls match a single stream()", {
  fit <- USAccDeaths_tbl %>%
    head(-12) %>%
    model(theta = THETA(value))
  streamed <- fit %>% stream(USAccDeaths_tbl %>% tail(12))

  # One observation at a time exercises the seasonal index rotation
  chained <- fit
  for (i in seq_len(12)) {
    chained <- chained %>% stream(USAccDeaths_tbl[NROW(USAccDeaths_tbl) - 12 + i, ])
  }
  expect_equal(fitted(chained)[[".fitted"]], fitted(streamed)[[".fitted"]])
  expect_equal(residuals(chained)[[".resid"]], residuals(streamed)[[".resid"]])
  expect_equal(glance(chained), glance(streamed))
  expect_equal(chained$theta[[1]]$fit$lT, streamed$theta[[1]]$fit$lT)
  expect_equal(chained$theta[[1]]$fit$season, streamed$theta[[1]]$fit$season)
  expect_equal(forecast(chained, h = 24)$.mean, forecast(streamed, h = 24)$.mean)

  # Uneven chunks (not a multiple of the seasonal period)
  chunked <- fit %>%
    stream(USAccDeaths_tbl %>% tail(12) %>% head(5)) %>%
    stream(USAccDeaths_tbl %>% tail(7))
  expect_equal(fitted(chunked)[[".fitted"]], fitted(streamed)[[".fitted"]])
  expect_equal(chunked$theta[[1]]$fit$season, streamed$theta[[1]]$fit$season)
  expect_equal(forecast(chunked, h = 24)$.mean, forecast(streamed, h = 24)$.mean)
})

test_that("stream.fable_theta handles missing values", {
  dat <- USAccDeaths_tbl
  dat$value[c(63, 68)] <- NA
  fit <- dat %>%
    head(-12) %>%
    model(theta = THETA(value))
  mdl <- fit$theta[[1]]$fit
  nw <- dat %>% tail(12)

  streamed <- fit %>% stream(nw)
  expected <- theta_manual_stream(mdl, nw$value)

  n_old <- NROW(dat) - 12
  expect_equal(fitted(streamed)[[".fitted"]][-seq_len(n_old)], expected$fitted)
  expect_equal(residuals(streamed)[[".resid"]][-seq_len(n_old)], expected$resid)
  expect_equal(which(is.na(residuals(streamed)[[".resid"]])), c(63, 68))
  expect_false(any(is.na(fitted(streamed)[[".fitted"]])))
  expect_equal(glance(streamed)$sigma2, expected$sigma2)
  expect_false(any(is.na(forecast(streamed, h = 6)$.mean)))
})

test_that("stream.fable_theta must start immediately after the trained data", {
  fit <- USAccDeaths_tbl %>%
    head(-12) %>%
    model(theta = THETA(value))
  expect_error(
    fit %>% stream(USAccDeaths_tbl %>% tail(11)),
    "must start one step beyond the end of"
  )
  expect_error(
    fit %>% stream(USAccDeaths_tbl %>% tail(13)),
    "must start one step beyond the end of"
  )
})
