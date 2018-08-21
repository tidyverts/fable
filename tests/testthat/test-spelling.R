context("test-spelling.R")

test_that("package spell check", {
  skip_on_cran()
  
  # Determine package source path by finding README.md
  package_dir <- dirname(list.files("../../", pattern = "README.md", recursive = TRUE, full.names = TRUE))
  skip_if(length(package_dir) != 1) 
  badspell <- spelling::spell_check_package(
    package_dir
  )
  
  expect_equal(NROW(badspell), 0, info = capture.output(print(badspell)))
})
