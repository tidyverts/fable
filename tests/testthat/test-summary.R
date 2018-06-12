context("test-summary.R")

test_that("Summarising mables / fables", {
  fit <- USAccDeaths %>% 
    as_tsibble() %>%
    ETS(value)
  
  fc <- fit %>% 
    forecast()
  
  fc_sum_default <- fc %>% summary
  
  expect_equal(dim(fc_sum), c(24, 4))
  expect_s3_class(fc_sum$`80%`, "hilo")
})
