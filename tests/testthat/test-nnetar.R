context("test-nnetar")

test_that("NNETAR", {
  airmiles <- as_tsibble(airmiles)
  air_fit <- airmiles %>% model(NNETAR(BoxCox(value, 0.15)))
  expect_output(print(air_fit), regexp = "NNAR\\(1,1\\)")
  
  # Test with xreg
  air_fit <- airmiles %>% model(NNETAR(BoxCox(value, 0.15) ~ trend()))
  air_fit <- airmiles %>% model(NNETAR(BoxCox(value, 0.15) ~ trend() + rnorm(length(index))))
  
  # Test simulations
  air_fit %>% 
    imitate(h=10, times = 100)
  
  # Test forecasts
  air_fit %>% 
    forecast(h=10, times = 100)
})
