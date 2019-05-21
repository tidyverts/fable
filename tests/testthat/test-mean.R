context("test-mean")

test_that("MEAN", {
  fit <- USAccDeaths_tbl %>% model(mean = MEAN(value))
  expect_identical(
    fitted(fit) %>% select(index, .fitted),
    USAccDeaths_tbl %>% transmute(.fitted = mean(value))
  )
  
  expect_identical(
    glance(fit)$sigma, sd(scale(USAccDeaths_tbl$value, scale = FALSE))
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
    fc$value, rep(mean(USAccDeaths_tbl$value), 3)
  )
})
