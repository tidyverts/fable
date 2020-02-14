context("test-nnetar")

airmiles <- as_tsibble(airmiles)
test_that("Automatic NNETAR selection", {
  air_fit <- airmiles %>% model(NNETAR(box_cox(value, 0.15)))
  expect_equal(model_sum(air_fit[[1]][[1]]), "NNAR(1,1)")

  # Test with xreg
  air_fit <- airmiles %>% model(NNETAR(box_cox(value, 0.15) ~ trend() + rnorm(length(index))))
  air_fit <- airmiles %>% model(NNETAR(box_cox(value, 0.15) ~ trend()))

  # Test simulations
  air_fit %>%
    generate(h = 10, times = 5)

  # Test forecasts
  fc_sim <- air_fit %>%
    forecast(h = 10, times = 5)
  fc_boot <- air_fit %>%
    forecast(h = 10, times = 5, bootstrap = TRUE)
  expect_equal(
    fc_sim$value,
    fc_boot$value,
    tolerance = 100
  )

  # Short series
  expect_output(
    UKLungDeaths[1:24, ] %>%
      model(NNETAR(mdeaths)) %>%
      report(),
    "NNAR\\(4,1,3\\)\\[12\\]"
  )
})

test_that("Manual NNETAR selection", {
  fit <- UKLungDeaths %>%
    model(NNETAR(mdeaths ~ AR(p = 3, P = 2)))
  expect_equal(model_sum(fit[[1]][[1]]), "NNAR(3,2,3)[12]")

  expect_equal(
    with(augment(fit), .fitted + .resid)[-(1:24)],
    UKLungDeaths$mdeaths[-(1:24)]
  )

  # Short series
  expect_warning(
    airmiles[1:5, ] %>%
      model(NNETAR(value ~ AR(10))),
    "Reducing number of lagged inputs due to short series"
  )
})


test_that("NNETAR with bad inputs", {
  expect_warning(
    airmiles[1:2, ] %>%
      model(NNETAR(value)),
    "Not enough data to fit a model"
  )

  expect_warning(
    airmiles %>%
      model(NNETAR(resp(rep_along(value, NA)))),
    "All observations are missing, a model cannot be estimated without data"
  )

  expect_warning(
    airmiles %>%
      model(NNETAR(resp(rep_along(value, 1)))),
    "Constant data, setting `AR\\(p=1, P=0\\)`, and `scale_inputs=FALSE`"
  )

  expect_warning(
    airmiles %>%
      model(NNETAR(value ~ rep_along(value, 1))),
    "Constant xreg column, setting `scale_inputs=FALSE`"
  )
})
