context("test-spelling.R")

test_that("package spell check", {
  badspell <- spelling::spell_check_package(file.path(system.file(package="fable"), ".."))
  expect_equal(NROW(badspell), 0, 
               info = capture.output(print(badspell)))
})
