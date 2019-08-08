test_that("barplotCounts ProteinExperiment", {

  PE <- testList[[1]]

  expect_is(barplotCounts(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'counts'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')
})

test_that("barplotCounts PeptideExperiment", {

  PE <- testList[[2]]

  expect_is(barplotCounts(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'counts'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')

})

test_that("barplotCounts ProteomicsExperiment", {

  PE <- testList[[3]]

  expect_is(barplotCounts(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotCounts(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'counts'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')
  expect_true('mode' %in% colnames(plotDf))
  expect_is(plotDf[,'mode'], 'factor')

})
