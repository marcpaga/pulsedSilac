context("test-coldata")

test_that("colData getter works", {

  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_silent(cd <- colData(testProtExp))
  expect_is(cd, 'DataFrame')
  expect_equal(nrow(cd), 4)
  expect_equal(ncol(cd), 3)

  expect_silent(cd <- colData(testPeptExp))
  expect_is(cd, 'DataFrame')
  expect_equal(nrow(cd), 4)
  expect_equal(ncol(cd), 3)

  expect_silent(cd <- colData(testPE))
  expect_is(cd, 'DataFrame')
  expect_equal(nrow(cd), 4)
  expect_equal(ncol(cd), 3)

})

test_that("colData setter works", {

  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_silent(colData(testProtExp) <- data.frame(a = 1:4))
  expect_silent(colData(testProtExp) <- data.frame(a = 1:4, b = 1:4))
  expect_error(colData(testProtExp) <- data.frame(a = 1:5))

  expect_silent(colData(testPeptExp) <- data.frame(a = 1:4))
  expect_silent(colData(testPeptExp) <- data.frame(a = 1:4, b = 1:4))
  expect_error(colData(testPeptExp) <- data.frame(a = 1:5))

  expect_silent(colData(testPE) <- data.frame(a = 1:4))
  expect_silent(colData(testPE) <- data.frame(a = 1:4, b = 1:4))
  expect_error(colData(testPE) <- data.frame(a = 1:5))


})
