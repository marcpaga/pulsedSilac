context("test-dim")

test_that("dim works", {

  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_equal(dim(testProtExp), c(3,4))

  expect_equal(dim(testPeptExp), c(5,4))

  expect_equivalent(dim(testPE), matrix(c(3, 4, 5, 4), ncol = 2, byrow = TRUE))


})
