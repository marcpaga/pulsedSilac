context("test-ncol")

test_that("ncol works", {

  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_equal(ncol(testProtExp), 4)
  expect_equal(ncol(testProtExp[1,]), 4)
  expect_equal(ncol(testProtExp[,1]), 1)

  expect_equal(ncol(testPeptExp), 4)
  expect_equal(ncol(testPeptExp[1,]), 4)
  expect_equal(ncol(testPeptExp[,1]), 1)

  expect_equal(ncol(testPE), 4)
  expect_equal(ncol(testPE[1,]), 4)
  expect_equal(ncol(testPE[,1]), 1)


})
