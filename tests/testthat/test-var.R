context("test-var")

test_that("Automatic VAR selection", {
  fit <- UKLungDeaths %>%
    model(VAR(vars(mdeaths, fdeaths) ~ AR(0:2)))

  expect_output(
    report(fit),
    "Model: VAR\\(2\\) w/ mean"
  )

  expect_equal(NROW(tidy(fit)), 10)

  expect_identical(
    dim(glance(fit)$sigma2[[1]]),
    rep(2L, 2)
  )

  expect_equal(
    matrix(with(augment(fit), .fitted + .resid), ncol = 2)[3:72, ],
    cbind(UKLungDeaths$mdeaths[3:72], UKLungDeaths$fdeaths[3:72])
  )
})

test_that("Univariate VAR", {
  fit <- UKLungDeaths %>%
    model(VAR(fdeaths ~ 0))
  expect_s3_class(
    forecast(fit),
    "fbl_ts"
  )
})

test_that("VAR with xregs", {
  fit <- UKLungDeaths %>%
    model(VAR(vars(mdeaths, fdeaths) ~ AR(2) + fourier(K = 5)))

  expect_equal(
    NROW(tidy(fit)),
    30
  )

  expect_s3_class(
    forecast(fit),
    "fbl_ts"
  )
})
