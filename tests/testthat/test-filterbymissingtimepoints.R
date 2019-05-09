context("test-filterMissingTimepoints")

test_that("filterByMissingTimepoints works", {

  PE <- testList[[1]]
  metaoptions(PE)[['conditionCol']] <- 'condition'
  assays(PE)[[2]][2, 1:2] <- NA
  assays(PE)[[2]][3, 3:4] <- NA

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 0,
                                               returnVector = TRUE))

  expect_equal(length(v), nrow(PE))
  expect_is(v, 'logical')
  expect_equal(v, c(T, F, F))

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'int',
                                               maxMissing = 0,
                                               returnVector = TRUE))
  expect_equal(v, c(T, T, T))

  assays(PE)[[2]][2, 2] <- 1
  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 1,
                                               returnVector = TRUE))
  expect_equal(v, c(T, T, F))

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 0,
                                               returnVector = TRUE,
                                               strict = TRUE))
  expect_equal(v, c(T, F, F))

})

test_that("overlapFeatures works", {

  PE <- testList[[1]]
  metaoptions(PE)[['conditionCol']] <- 'condition'
  assays(PE)[[2]][2, 1:2] <- NA
  assays(PE)[[2]][3, 3:4] <- NA

  expect_silent(v <- overlapFeatures(PE,
                                     assayName = 'ratio',
                                     maxMissing = 0,
                                     returnList = TRUE))

  expect_is(v, 'list')
  v <- unname(v)
  expect_equal(v, list(c(1, 3),
                       c(1, 2)))

  expect_silent(v <- overlapFeatures(PE,
                                     assayName = 'ratio',
                                     maxMissing = 2,
                                     returnList = TRUE))
  v <- unname(v)
  expect_equal(v, list(c(1, 2, 3),
                       c(1, 2, 3)))

  assays(PE)[[2]][2, 2] <- 1
  expect_silent(v <- overlapFeatures(PE,
                                     assayName = 'ratio',
                                     maxMissing = 1,
                                     returnList = TRUE))
  v <- unname(v)
  expect_equal(v, list(c(1, 2, 3),
                       c(1, 2)))

})
