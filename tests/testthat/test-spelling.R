context("test-spelling.R")

test_that("package spell check", {
  badspell <- spelling::spell_check_package(
    # Determine package source path by finding README.md
    dirname(list.files("../../", pattern = "README.md", recursive = TRUE, full.names = TRUE))
  )
  expect_equal(NROW(badspell), 0, 
               info = capture.output(print(badspell)))
})
