context("test-mean")

test_that("MEAN", {
  fable_fit <- USAccDeaths_tbl %>% model(mean = MEAN(value))
  expect_identical(
    fitted(fable_fit) %>% select(index, .fitted),
    USAccDeaths_tbl %>% transmute(.fitted = mean(value))
  )
})
