context("test-multivariate.R")

test_that("multiple univariate", {
  fit <- tsibbledata::UKLungDeaths %>%
    gather(type, deaths, -index) %>%
    ETS(deaths)
  
  expect_equal(fit$type, c("fdeaths", "mdeaths"))
  
  fc <- fit %>% forecast
  
  expect_equal(fc$type, c("fdeaths", "mdeaths"))
})
