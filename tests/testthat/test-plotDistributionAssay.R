test_that("plotDistributionAssay works", {

  PE <- testList[[1]]
  expect_error(plotDistributionAssay(x = PE,
                                     assayName = 'asdf'))

  expect_error(plotDistributionAssay(x = PE,
                                     assayName = 'ratio',
                                     plotType = 'asdf'))

  expect_silent(p <- plotDistributionAssay(x = PE,
                                     assayName = 'ratio',
                                     plotType = 'density'))
  expect_is(p, 'ggplot')
  expect_equal(sapply(p$layers, function(x) class(x$geom)[1]),
               'GeomDensityRidges')

  expect_silent(p <- plotDistributionAssay(x = PE,
                                           assayName = 'ratio',
                                           plotType = 'boxplot'))
  expect_is(p, 'ggplot')
  expect_equal(sapply(p$layers, function(x) class(x$geom)[1]), 'GeomBoxplot')

  expect_silent(p <- plotDistributionAssay(x = PE,
                                           assayName = 'ratio',
                                           plotType = 'boxplot',
                                           returnDataFrame = TRUE))
  expect_is(p, 'data.frame')
  expect_named(p, c('ratio', 'time', 'condition'))
  expect_equal(nrow(p), length(assays(PE)[[1]]))


  PE <- testList[[3]]
  expect_silent(p <- plotDistributionAssay(x = PE,
                                           assayName = 'ratio',
                                           plotType = 'boxplot',
                                           mode = 'protein',
                                           returnDataFrame = TRUE))
  expect_is(p, 'data.frame')
  expect_named(p, c('ratio', 'time', 'condition'))
  expect_equal(nrow(p), length(assaysProt(PE)[[1]]))

  expect_silent(p <- plotDistributionAssay(x = PE,
                                           assayName = 'ratio',
                                           plotType = 'boxplot',
                                           mode = 'peptide',
                                           returnDataFrame = TRUE))
  expect_is(p, 'data.frame')
  expect_named(p, c('ratio', 'time', 'condition'))
  expect_equal(nrow(p), length(assaysPept(PE)[[1]]))


})
