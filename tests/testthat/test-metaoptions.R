context("test-metaoptions")

test_that("metaoptions ProteinExperiment", {

  ## test objects
  testProtExp <- testList[[1]]

  ## getter
  expect_silent(mo <- metaoptions(testProtExp))
  expect_equal(length(mo), 4)
  expect_equal(names(mo), c('conditionCol',
                            'timeCol',
                            'replicateIntCol',
                            'replicateTimeCol'))

  ## setter
  expect_silent(metaoptions(testProtExp)[['conditionCol']] <- 'a')
  expect_silent(metaoptions(testProtExp)[['timeCol']] <- 3)
  expect_silent(metaoptions(testProtExp)[['timeCol']] <- NA)
  expect_error(metaoptions(testProtExp)[['timeCol']] <- TRUE)
  expect_error(metaoptions(testProtExp)[['timeCol']] <- as.factor('a'))

})

test_that("metaoptions PeptideExperiment", {

  ## test objects
  testPeptExp <- testList[[2]]

  ## getter
  expect_silent(mo <- metaoptions(testPeptExp))
  expect_equal(length(mo), 5)
  expect_equal(names(mo), c('conditionCol',
                            'timeCol',
                            'replicateIntCol',
                            'replicateTimeCol',
                            'proteinCol'))

  ## setter
  expect_silent(metaoptions(testPeptExp)[['conditionCol']] <- 'a')
  expect_silent(metaoptions(testPeptExp)[['timeCol']] <- 3)
  expect_silent(metaoptions(testPeptExp)[['timeCol']] <- NA)
  expect_error(metaoptions(testPeptExp)[['timeCol']] <- TRUE)
  expect_error(metaoptions(testPeptExp)[['timeCol']] <- as.factor('a'))

})


test_that("metaoptions ProteomicsExperiment", {

  ## test objects
  testPE <- testList[[3]]

  ## getter
  expect_silent(mo <- metaoptions(testPE))
  expect_equal(length(mo), 9)
  expect_true(all(names(mo) %in% c('conditionCol',
                                   'timeCol',
                                   'replicateIntCol',
                                   'replicateTimeCol',
                                   'idColProt',
                                   'idColPep',
                                   'linkedSubset',
                                   'subsetMode',
                                   'proteinCol')))

  ## setter
  expect_silent(metaoptions(testPE)[['conditionCol']] <- 'a')
  expect_silent(metaoptions(testPE)[['timeCol']] <- 3)
  expect_silent(metaoptions(testPE)[['timeCol']] <- NA)
  expect_error(metaoptions(testPE)[['timeCol']] <- TRUE)
  expect_error(metaoptions(testPE)[['timeCol']] <- as.factor('a'))
  expect_silent(metaoptions(testPE)[['linkedSubset']] <- FALSE)
  expect_error(metaoptions(testPE)[['linkedSubset']] <- 3)
  expect_error(metaoptions(testPE)[['linkedSubset']] <- '3')
  expect_error(metaoptions(testPE)[['subsetMode']] <- FALSE)
  expect_error(metaoptions(testPE)[['subsetMode']] <- 3)
  expect_error(metaoptions(testPE)[['subsetMode']] <- '3')
  expect_silent(metaoptions(testPE)[['subsetMode']] <- 'peptide')
})
