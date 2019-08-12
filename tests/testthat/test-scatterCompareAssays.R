test_that("scatterCompareAssays works", {

  PE <- testList[[1]]
  metaoptions(PE)[['conditionCol']] <- 'condition'
  metaoptions(PE)[['timeCol']] <- 'time'

  expect_error(scatterCompareAssays(PE,
                                    conditions = c('1', '3'),
                                    assayName = 'int'))

  expect_error(scatterCompareAssays(PE,
                                    conditions = c('1', '2'),
                                    assayName = 'asdf'))

  expect_silent(p <- scatterCompareAssays(PE,
                                          conditions = c('1', '2'),
                                          assayName = 'int'))
  expect_is(p, 'ggplot')

  expect_silent(p <- scatterCompareAssays(PE,
                                          conditions = c('1', '2'),
                                          assayName = 'int',
                                          returnDataFrame = TRUE))
  expect_is(p, 'data.frame')
  expect_equal(colnames(p)[3], 'Time')
  expect_equal(colnames(p)[1:2],  c('1', '2'))
  expect_equal(nrow(p), length(assays(PE)[[1]]))

  PE <- testList[[3]]
  metaoptions(PE)[['conditionCol']] <- 'condition'
  metaoptions(PE)[['timeCol']] <- 'time'

  expect_silent(p <- scatterCompareAssays(PE,
                                          conditions = c('1', '2'),
                                          assayName = 'int',
                                          mode = 'protein',
                                          returnDataFrame = TRUE))

  expect_is(p, 'data.frame')
  expect_equal(colnames(p)[3], 'Time')
  expect_equal(colnames(p)[1:2],  c('1', '2'))
  expect_equal(nrow(p), length(assaysProt(PE)[[1]]))

  expect_silent(p <- scatterCompareAssays(PE,
                                          conditions = c('1', '2'),
                                          assayName = 'int',
                                          mode = 'peptide',
                                          returnDataFrame = TRUE))

  expect_is(p, 'data.frame')
  expect_equal(colnames(p)[3], 'Time')
  expect_equal(colnames(p)[1:2],  c('1', '2'))
  expect_equal(nrow(p), length(assaysPept(PE)[[1]]))


})
