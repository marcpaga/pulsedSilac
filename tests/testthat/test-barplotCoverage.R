test_that("barplotTimeCoverage ProteinExperiment", {

  PE <- testList[[1]]

  expect_is(barplotTimeCoverage(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'Freq'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')
})

test_that("barplotTimeCoverage PeptideExperiment", {

  PE <- testList[[2]]

  expect_is(barplotTimeCoverage(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'Freq'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')

})

test_that("barplotTimeCoverage ProteomicsExperiment", {

  PE <- testList[[3]]

  expect_is(barplotTimeCoverage(x = PE, assayName = 'int'), 'ggplot')
  expect_is(barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE), 'data.frame')

  plotDf <- barplotTimeCoverage(x = PE, assayName = 'int',
                          returnDataFrame = TRUE)

  expect_true('counts' %in% colnames(plotDf))
  expect_is(plotDf[,'Freq'], 'integer')
  expect_true('condition' %in% colnames(plotDf))
  expect_is(plotDf[,'condition'], 'factor')
  expect_true('mode' %in% colnames(plotDf))
  expect_is(plotDf[,'mode'], 'factor')

})
