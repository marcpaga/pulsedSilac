context("test-assayNames")

test_that("assayNames getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(assayNames(testProtExp), 'character')
  expect_silent(an <- assayNames(testProtExp))
  expect_equal(length(an), 2)

  ## PeptideExperiment
  expect_is(assayNames(testPeptExp), 'character')
  expect_silent(an <- assays(testPeptExp))
  expect_equal(length(an), 2)

  ## ProteomicsExperiment
  expect_is(assayNames(testPE), 'list')
  expect_silent(an <- assays(testPE))
  expect_equal(length(an), 2)
  expect_equal(names(an), c('protein', 'peptide'))
  expect_equal(names(an[[1]]), c('int', 'ratio'))
  expect_equal(names(an[[2]]), c('int', 'ratio'))

})

test_that("assayNames setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment

  expect_silent(assayNames(testProtExp)[[1]] <- 'int2')
  expect_equal(assayNames(testProtExp), c('int2', 'ratio'))
  expect_error(assayNames(testProtExp)[[3]] <- 'int')

  ## Peptide Experiment
  expect_silent(assayNames(testPeptExp)[[1]] <- 'int2')
  expect_equal(assayNames(testPeptExp), c('int2', 'ratio'))
  expect_error(assayNames(testPeptExp)[[3]] <- 'int')


})


test_that("assayNamesProt getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_silent(assayNamesProt(testProtExp)[[1]] <- 'int2')
  expect_equal(assayNamesProt(testProtExp), c('int2', 'ratio'))
  expect_error(assayNamesProt(testProtExp)[[3]] <- 'int')

  ## ProteomicsExperiment
  expect_silent(assayNamesProt(testPE)[[1]] <- 'int2')
  expect_equal(assayNamesProt((testPE)), c('int2', 'ratio'))
  expect_error(assayNamesProt((testPE))[[3]] <- 'int')

})


test_that("assayNamesProt setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_silent(assayNamesProt(testProtExp)[[1]] <- 'int2')
  expect_equal(assayNamesProt(testProtExp), c('int2', 'ratio'))
  expect_error(assayNamesProt(testProtExp)[[3]] <- 'int')


  ## ProteomicsExperiment
  expect_silent(assayNamesProt(testPE)[[1]] <- 'int2')
  expect_equal(assayNamesProt(testPE), c('int2', 'ratio'))
  expect_error(assayNamesProt(testPE)[[3]] <- 'int')



})


test_that("assayNamesPept getter", {

  ## test objects
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## PeptideExperiment
  expect_silent(assayNamesPept(testPeptExp)[[1]] <- 'int2')
  expect_equal(assayNamesPept(testPeptExp), c('int2', 'ratio'))
  expect_error(assayNamesPept(testPeptExp)[[3]] <- 'int')

  ## ProteomicsExperiment
  expect_silent(assayNamesPept(testPE)[[1]] <- 'int2')
  expect_equal(assayNamesPept((testPE)), c('int2', 'ratio'))
  expect_error(assayNamesPept((testPE))[[3]] <- 'int')

})


test_that("assayNamesPept setter", {

  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## PeptideExperiment
  expect_silent(assayNamesPept(testPeptExp)[[1]] <- 'int2')
  expect_equal(assayNamesPept(testPeptExp), c('int2', 'ratio'))
  expect_error(assayNamesPept(testPeptExp)[[3]] <- 'int')


  ## ProteomicsExperiment
  expect_silent(assayNamesPept(testPE)[[1]] <- 'int2')
  expect_equal(assayNamesPept(testPE), c('int2', 'ratio'))
  expect_error(assayNamesPept(testPE)[[3]] <- 'int')
})

