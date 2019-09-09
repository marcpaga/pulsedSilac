context("test-linkerdf")


test_that("linkerDf getter", {

  testPE <- testList[[3]]

  expect_silent(df <- linkerDf(testPE))
  expect_equal(nrow(df), 5)
  expect_equal(ncol(df), 4)

})

test_that("linkerDf setter", {

  testPE <- testList[[3]]

  df <- linkerDf(testPE)
  expect_silent(linkerDf(testPE) <- df)
  expect_error(linkerDf(testPE) <- df[,1:3])
  expect_error(linkerDf(testPE) <- df[1:2,])
  metadata(testPE)[['linkedSubset']] <- FALSE
  expect_silent(linkerDf(testPE) <- df[1:2,])

})
