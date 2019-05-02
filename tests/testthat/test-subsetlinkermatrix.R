context("test-subsetlinkerdf")

test_that("subsetlinkerdf works", {

  protIDs <- c(1:4)
  pepIDs <- c(1:6)
  protToPep <- list(c(1, 2), c(3, 4), c(5, 6), c(5, 6))
  pepToProt <- list(1, 1, 2, 2, c(3, 4), c(3, 4))

  linkerMat <- buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt)

  expect_error(subsetLinkerDf(protIDs = c(1:2)))
  expect_silent(subsetLinkerDf(x = linkerMat,
                                   protIDs = c(1:2)))
  expect_error(subsetLinkerDf(x = linkerMat,
                                   protIDs = c(1:2),
                                  protRows = c(1:4)))
  expect_error(subsetLinkerDf(x = linkerMat,
                                  protIDs = c(1:2),
                                  protRows = c(1:4),
                                  pepRows = c(1)))
  expect_error(subsetLinkerDf(x = linkerMat,
                                  protIDs = c(1:2),
                                  protRows = c(1:4),
                                  pepRows = c(1),
                                  pepIDs = c(1)))

  slm <- subsetLinkerDf(x = linkerMat, protIDs = c(1:2))
  expect_equal(nrow(slm), 4)
  expect_equal(ncol(slm), 4)
  slm <- subsetLinkerDf(x = linkerMat, protIDs = 3)
  expect_equal(nrow(slm), 4)
  slm <- subsetLinkerDf(x = linkerMat, protRows = c(2:3))
  expect_equal(nrow(slm), 6)
  slm <- subsetLinkerDf(x = linkerMat, pepIDs = c(4, 5, 6))
  expect_equal(nrow(slm), 6)
  slm <- subsetLinkerDf(x = linkerMat, pepRows = 5)
  expect_equal(nrow(slm), 4)

})
