context("test-nrow")

test_that("nrow works", {

  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_equal(nrow(testProtExp), 3)
  expect_equal(nrow(testProtExp[1,]), 1)

  expect_equal(nrow(testPeptExp), 5)
  expect_equal(nrow(testPeptExp[1,]), 1)

  expect_equivalent(nrow(testPE), c(3, 5))
  expect_equivalent(nrow(testPE[1,]), c(1, 2))
})
