context("test-rowdata")


test_that("rowData getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(rowData(testProtExp), 'DataFrame')
  expect_silent(rd <- rowData(testProtExp))
  expect_equal(nrow(rd), nrow(testProtExp))
  expect_equal(ncol(rd), 2)

  ## PeptideExperiment
  expect_is(rowData(testPeptExp), 'DataFrame')
  expect_silent(rd <- rowData(testPeptExp))
  expect_equal(nrow(rd), nrow(testPeptExp))
  expect_equal(ncol(rd), 2)

  ## ProteomicsExperiment
  expect_is(rowData(testPE), 'list')
  expect_silent(rd <- rowData(testPE))
  expect_equal(length(rd), 2)
  expect_equal(names(rd), c('protein', 'peptide'))

})

test_that("rowData setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  new.rd <- data.frame(id = LETTERS[4:6])
  expect_silent(rowData(testProtExp) <- new.rd)
  expect_equal(ncol(rowData(testProtExp)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowData(testProtExp) <- new.rd)
  new.rd <- data.frame(id = LETTERS[4:10])
  expect_error(rowData(testProtExp) <- new.rd)

  ## Peptide Experiment
  new.rd <- data.frame(id = LETTERS[6:10])
  expect_silent(rowData(testPeptExp) <- new.rd)
  expect_equal(ncol(rowData(testPeptExp)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowData(testPeptExp) <- new.rd)
  new.rd <- data.frame(id = LETTERS[1:10])
  expect_error(rowData(testPeptExp) <- new.rd)

  ## ProteomicsExperiment
  new.rd <- list(data.frame(id = LETTERS[4:6]),
                 data.frame(id = LETTERS[6:10]))
  expect_silent(rowData(testPE) <- new.rd)
  new.rd <- list(data.frame(id = LETTERS[2:6]),
                 data.frame(id = LETTERS[6:10]))
  expect_error(rowData(testPE) <- new.rd)
  new.rd <- list(data.frame(id = LETTERS[4:6]),
                 data.frame(id = LETTERS[8:10]))
  expect_error(rowData(testPE) <- new.rd)

})


test_that("rowDataProt getter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_is(rowDataProt(testProtExp), 'DataFrame')
  expect_silent(rd <- rowDataProt(testProtExp))
  expect_equal(nrow(rd), nrow(testProtExp))
  expect_equal(ncol(rd), 2)

  ## ProteomicsExperiment
  expect_is(rowDataProt(testPE), 'DataFrame')
  expect_silent(rd <- rowDataProt(testPE))
  expect_equal(nrow(rd), unname(nrow(testPE)[1]))
  expect_equal(ncol(rd), 2)

})


test_that("rowDataProt setter", {

  ## test objects
  testProtExp <- testList[[1]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  new.rd <- data.frame(id = LETTERS[4:6])
  expect_silent(rowDataProt(testProtExp) <- new.rd)
  expect_equal(ncol(rowData(testProtExp)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowDataProt(testProtExp) <- new.rd)
  new.rd <- data.frame(id = LETTERS[4:10])
  expect_error(rowDataProt(testProtExp) <- new.rd)

  ## ProteomicsExperiment
  new.rd <- data.frame(id = LETTERS[4:6])
  expect_silent(rowDataProt(testPE) <- new.rd)
  expect_equal(ncol(rowDataProt(testPE)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowDataProt(testPE) <- new.rd)
  new.rd <- data.frame(id = LETTERS[4:10])
  expect_error(rowDataProt(testPE) <- new.rd)

})


test_that("rowDataPept getter", {

  ## test objects
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## PeptideExperiment
  expect_is(rowDataPept(testPeptExp), 'DataFrame')
  expect_silent(rd <- rowDataPept(testPeptExp))
  expect_equal(nrow(rd), nrow(testPeptExp))
  expect_equal(ncol(rd), 2)

  ## ProteomicsExperiment
  expect_is(rowDataPept(testPE), 'DataFrame')
  expect_silent(rd <- rowDataPept(testPE))
  expect_equal(nrow(rd), unname(nrow(testPE)[2]))
  expect_equal(ncol(rd), 2)

})


test_that("rowDataPept setter", {

  ## test objects
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## Peptide Experiment
  new.rd <- data.frame(id = LETTERS[6:10])
  expect_silent(rowDataPept(testPeptExp) <- new.rd)
  expect_equal(ncol(rowDataPept(testPeptExp)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowDataPept(testPeptExp) <- new.rd)
  new.rd <- data.frame(id = LETTERS[1:10])
  expect_error(rowDataPept(testPeptExp) <- new.rd)

  ## ProteomicsExperiment
  new.rd <- data.frame(id = LETTERS[6:10])
  expect_silent(rowDataPept(testPE) <- new.rd)
  expect_equal(ncol(rowDataPept(testPE)), 1)

  new.rd <- data.frame(id = LETTERS[4:5])
  expect_error(rowDataPept(testPE) <- new.rd)
  new.rd <- data.frame(id = LETTERS[1:10])
  expect_error(rowDataPept(testPE) <- new.rd)

})
