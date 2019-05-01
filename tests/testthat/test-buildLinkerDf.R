context("test-buildLinkerMatrix")

test_that("function works redundant parameters", {

  protIDs <- c(1:3)
  pepIDs <- c(1:6)
  protToPep <- list(c(1,2), c(3,4), c(5,6))
  pepToProt <- list(1, 1, 2, 2, 3, 3)

  linkerMat <- buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt)

  expect_equal(nrow(linkerMat), 6)
  expect_equal(ncol(linkerMat), 4)
  expect_equal(sum(is.na(linkerMat)), 0)

  expect_is(linkerMat[,1], c('numeric', 'character'))
  expect_is(linkerMat[,2], c('numeric', 'character'))
  expect_is(linkerMat[,3], 'numeric')
  expect_is(linkerMat[,4], 'numeric')

})

test_that("function works reduced params 1", {

  protIDs <- c(1:3)
  pepIDs <- c(1:6)
  protToPep <- list(c(1,2), c(3,4), c(5,6))

  linkerMat <- buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep =  protToPep)

  expect_equal(nrow(linkerMat), 6)
  expect_equal(ncol(linkerMat), 4)
  expect_equal(sum(is.na(linkerMat)), 0)

  expect_is(linkerMat[,1], c('numeric', 'character'))
  expect_is(linkerMat[,2], c('numeric', 'character'))
  expect_is(linkerMat[,3], 'numeric')
  expect_is(linkerMat[,4], 'numeric')

})

test_that("function works reduced params 2", {

  protIDs <- c(1:3)
  pepIDs <- c(1:6)
  pepToProt <- list(1, 1, 2, 2, 3, 3)

  linkerMat <- buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             pepToProt = pepToProt)

  expect_equal(nrow(linkerMat), 6)
  expect_equal(ncol(linkerMat), 4)
  expect_equal(sum(is.na(linkerMat)), 0)

  expect_is(linkerMat[,1], c('numeric', 'character'))
  expect_is(linkerMat[,2], c('numeric', 'character'))
  expect_is(linkerMat[,3], 'numeric')
  expect_is(linkerMat[,4], 'numeric')

})


test_that("function raises proper errors", {

  protIDs <- c(1:2)
  pepIDs <- c(1:6)
  protToPep <- list(c(1,2), c(3,4), c(5,6))
  pepToProt <- list(1, 1, 2, 2, 3, 3)

  expect_error(buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'The ID vectors lengths do not match with the lists lengths')

  protIDs <- c(1:3)
  pepIDs <- c(1:5)
  protToPep <- list(c(1,2), c(3,4), c(5,6))
  pepToProt <- list(1, 1, 2, 2, 3, 3)

  expect_error(buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'The ID vectors lengths do not match with the lists lengths')


  expect_error(buildLinkerDf(protIDs = protIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'Both protIDs and pepIDs are necessary')

  expect_error(buildLinkerDf(pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'Both protIDs and pepIDs are necessary')

  protIDs <- c(1,2,3,3)
  pepIDs <- c(1:6)

  expect_error(buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'Not all ids are unique')

  protIDs <- c(1,2,3,3)
  pepIDs <- c(1:6, 6)

  expect_error(buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'Not all ids are unique')

  protIDs <- c(1,2,3)
  pepIDs <- c(1:6)
  protToPep <- list(c(1,2), c(3,4), c(5,6))
  pepToProt <- list(1, 1, 2, 2, 3, 4)

  expect_error(buildLinkerDf(protIDs = protIDs,
                             pepIDs = pepIDs,
                             protToPep = protToPep,
                             pepToProt = pepToProt),
               'There are some inconsistencies with the data given')
})

