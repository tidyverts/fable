context("test-multivariate.R")

test_that("multiple univariate", {
  fit <- tsibbledata::UKLungDeaths %>%
    gather(type, deaths, -index) %>%
    ETS(deaths)
  
  expect_equal(fit$type, c("fdeaths", "mdeaths"))
  expect_s3_class(fit$model, "lst_mdl")
  fc <- fit %>% forecast
  
  expect_equal(fc$type, c("fdeaths", "mdeaths"))
  expect_s3_class(fc$forecast, "lst_fc")
})
