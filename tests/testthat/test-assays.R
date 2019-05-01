context("test-assays")

test_that("assays getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(assays(testProtExp), 'SimpleList')
  expect_silent(as <- assays(testProtExp))
  expect_equal(length(as), 2)

  ## PeptideExperiment
  expect_is(assays(testPeptExp), 'SimpleList')
  expect_silent(as <- assays(testPeptExp))
  expect_equal(length(as), 2)

  ## ProteomicsExperiment
  expect_is(assays(testPE), 'list')
  expect_silent(as <- assays(testPE))
  expect_equal(length(as), 2)
  expect_equal(names(as), c('protein', 'peptide'))
  expect_equal(names(as[[1]]), c('int', 'ratio'))
  expect_equal(names(as[[2]]), c('int', 'ratio'))

})

test_that("assays setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_silent(assays(testProtExp)[[3]] <- new.as)
  expect_equal(length(assays(testProtExp)), 3)
  new.as <- matrix(data = 12:1, ncol = 4)
  expect_silent(assays(testProtExp)[[3]] <- new.as)
  expect_equal(length(assays(testProtExp)), 3)
  expect_equivalent(assays(testProtExp)[[3]], new.as)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assays(testProtExp)[[3]] <- new.as)
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_error(assays(testProtExp)[[5]] <- new.as)

  ## Peptide Experiment
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_silent(assays(testPeptExp)[[3]] <- new.as)
  expect_equal(length(assays(testPeptExp)), 3)
  new.as <- matrix(data = 20:1, ncol = 4)
  expect_silent(assays(testPeptExp)[[3]] <- new.as)
  expect_equal(length(assays(testPeptExp)), 3)
  expect_equivalent(assays(testPeptExp)[[3]], new.as)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assays(testPeptExp)[[3]] <- new.as)
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_error(assays(testPeptExp)[[5]] <- new.as)

})


test_that("assaysProt getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(assaysProt(testProtExp), 'SimpleList')
  expect_silent(as <- assaysProt(testProtExp))
  expect_equal(length(as), 2)

  ## ProteomicsExperiment
  expect_is(assaysProt(testPE), 'SimpleList')
  expect_silent(as <- assaysProt(testPE))
  expect_equal(length(as), 2)

})


test_that("assaysProt setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_silent(assaysProt(testProtExp)[[3]] <- new.as)
  expect_equal(length(assaysProt(testProtExp)), 3)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assaysProt(testProtExp)[[3]] <- new.as)
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_error(assaysProt(testProtExp)[[5]] <- new.as)

  ## ProteomicsExperiment
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_silent(assaysProt(testPE)[[3]] <- new.as)
  expect_equal(length(assaysProt(testPE)), 3)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assaysProt(testPE)[[3]] <- new.as)
  new.as <- matrix(data = 1:12, ncol = 4)
  expect_error(assaysProt(testPE)[[5]] <- new.as)


})


test_that("assaysPept getter", {

  ## test objects
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(assaysPept(testPeptExp), 'SimpleList')
  expect_silent(as <- assaysPept(testPeptExp))
  expect_equal(length(as), 2)

  ## ProteomicsExperiment
  expect_is(assaysPept(testPE), 'SimpleList')
  expect_silent(as <- assaysPept(testPE))
  expect_equal(length(as), 2)

})


test_that("assaysPept setter", {

  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_silent(assaysPept(testPeptExp)[[3]] <- new.as)
  expect_equal(length(assaysPept(testPeptExp)), 3)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assaysPept(testPeptExp)[[3]] <- new.as)
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_error(assaysPept(testPeptExp)[[5]] <- new.as)

  ## ProteomicsExperiment
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_silent(assaysPept(testPE)[[3]] <- new.as)
  expect_equal(length(assaysPept(testPE)), 3)

  new.as <- matrix(data = 1:12, ncol = 3)
  expect_error(assaysPept(testPE)[[3]] <- new.as)
  new.as <- matrix(data = 1:20, ncol = 4)
  expect_error(assaysPept(testPE)[[5]] <- new.as)

})

