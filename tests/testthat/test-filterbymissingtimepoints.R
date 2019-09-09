context("test-filterMissingTimepoints")

test_that("filterByMissingTimepoints works", {

  PE <- testList[[1]]
  metadata(PE)[['conditionCol']] <- 'condition'
  assays(PE)[[2]][2, 1:2] <- NA
  assays(PE)[[2]][3, 3:4] <- NA

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 0,
                                               returnVector = TRUE))

  expect_equal(length(v), nrow(PE))
  expect_is(v, 'logical')
  expect_equal(v, c(TRUE, FALSE, FALSE))

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'int',
                                               maxMissing = 0,
                                               returnVector = TRUE))
  expect_equal(v, c(TRUE, TRUE, TRUE))

  assays(PE)[[2]][2, 2] <- 1
  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 1,
                                               returnVector = TRUE))
  expect_equal(v, c(TRUE, TRUE, FALSE))

  expect_silent(v <- filterByMissingTimepoints(PE,
                                               assayName = 'ratio',
                                               maxMissing = 0,
                                               returnVector = TRUE,
                                               strict = TRUE))
  expect_equal(v, c(TRUE, FALSE, FALSE))

})

test_that("upsetTimeCoverage works", {

  PE <- testList[[1]]
  metadata(PE)[['conditionCol']] <- 'condition'
  assays(PE)[[2]][2, 1:2] <- NA
  assays(PE)[[2]][3, 3:4] <- NA

  expect_silent(v <- upsetTimeCoverage(PE,
                                     assayName = 'ratio',
                                     maxMissing = 0,
                                     returnList = TRUE))

  expect_is(v, 'list')
  v <- unname(v)
  expect_equal(v, list(c(1, 3),
                       c(1, 2)))

  expect_silent(v <- upsetTimeCoverage(PE,
                                     assayName = 'ratio',
                                     maxMissing = 2,
                                     returnList = TRUE))
  v <- unname(v)
  expect_equal(v, list(c(1, 2, 3),
                       c(1, 2, 3)))

  assays(PE)[[2]][2, 2] <- 1
  expect_silent(v <- upsetTimeCoverage(PE,
                                     assayName = 'ratio',
                                     maxMissing = 1,
                                     returnList = TRUE))
  v <- unname(v)
  expect_equal(v, list(c(1, 2, 3),
                       c(1, 2)))

})
