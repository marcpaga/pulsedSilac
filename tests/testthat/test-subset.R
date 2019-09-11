context("test-subset")

test_that("subset operator works", {

  protmatrix <- list(matrix(1:9, nrow = 3, ncol = 3))
  pepmatrix <- list(matrix(1:12, nrow = 4, ncol = 3))
  protdata <- data.frame(a = 1:3)
  pepdata <- data.frame(b = 1:4)
  linkM <- buildLinkerDf(protIDs = c(1:3),
                         pepIDs = c(1:4),
                         protToPep = list('1' = c(1),
                                          '2' = c(2),
                                          '3' = c(3, 4)))

  coldata <- data.frame(experiment = c('A', 'B', 'C'))

  protExp <- SilacProteinExperiment(assays = protmatrix,
                               rowData = protdata,
                               colData = coldata)
  peptExp <- SilacPeptideExperiment(assays = pepmatrix,
                               rowData = pepdata,
                               colData = coldata)


  PE <- SilacProteomicsExperiment(SilacProteinExperiment = protExp,
                                  SilacPeptideExperiment = peptExp,
                             colData = coldata,
                             linkerDf = linkM,
                             idColProt = 'a',
                             idColPept = 'b')

  metadata(PE)[['a']] <- 1

  ## protein level
  # test the rows
  expect_equal(as.vector(dim(PE[1, ])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[2, ])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[3, ])), c(1, 2, 3, 3))
  expect_error(PE[4, ])
  expect_equal(as.vector(dim(PE[1:3, ])), c(3, 4, 3, 3))
  expect_equal(as.vector(dim(PE[2:3, ])), c(2, 3, 3, 3))
  expect_error(PE[2:5, ])

  # test the columns
  expect_equal(as.vector(dim(PE[, 1])), c(3, 4, 1, 1))
  expect_equal(as.vector(dim(PE[, 2])), c(3, 4, 1, 1))
  expect_equal(as.vector(dim(PE[, 3])), c(3, 4, 1, 1))
  expect_error(PE[, 4])
  expect_equal(as.vector(dim(PE[, 1:3])), c(3, 4, 3, 3))
  expect_equal(as.vector(dim(PE[, 2:3])), c(3, 4, 2, 2))
  expect_error(PE[, 2:5])

  # test both
  expect_equal(as.vector(dim(PE[1,1])), c(1, 1, 1, 1))
  expect_error(PE[5, 1])
  expect_equal(as.vector(dim(PE[1:3,1])), c(3, 4, 1, 1))
  expect_equal(as.vector(dim(PE[1:3,1:3])), c(3, 4, 3, 3))

  ## peptide level
  metaoptions(PE)[['subsetMode']] <- 'peptide'
  expect_equal(as.vector(dim(PE[1, ])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[2, ])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[3, ])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[4, ])), c(1, 1, 3, 3))
  expect_error(PE[5, ])
  expect_equal(as.vector(dim(PE[1:3, ])), c(3, 3, 3, 3))
  expect_equal(as.vector(dim(PE[1:4, ])), c(3, 4, 3, 3))
  expect_equal(as.vector(dim(PE[2:4, ])), c(2, 3, 3, 3))
  expect_equal(as.vector(dim(PE[2:3, ])), c(2, 2, 3, 3))
  expect_error(PE[2:5, ])


  # test the columns
  expect_equal(as.vector(dim(PE[, 1])), c(3, 4, 1, 1))
  expect_equal(as.vector(dim(PE[, 2])), c(3, 4, 1, 1))
  expect_equal(as.vector(dim(PE[, 3])), c(3, 4, 1, 1))
  expect_error(PE[, 4])
  expect_equal(as.vector(dim(PE[, 1:3])), c(3, 4, 3, 3))
  expect_equal(as.vector(dim(PE[, 2:3])), c(3, 4, 2, 2))
  expect_error(PE[, 2:5])

  # test both
  expect_equal(as.vector(dim(PE[1,1])), c(1, 1, 1, 1))
  expect_error(PE[5, 1])
  expect_equal(as.vector(dim(PE[1:3,1])), c(3, 3, 1, 1))
  expect_equal(as.vector(dim(PE[1:3,1:3])), c(3, 3, 3, 3))
  expect_equal(as.vector(dim(PE[2,1:3])), c(1, 1, 3, 3))
  expect_equal(as.vector(dim(PE[2:3,1:3])), c(2, 2, 3, 3))

})


test_that("subset method works", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  ## ProteinExperiment
  expect_equal(dim(subset(testProtExp, id == 'A')), c(1, 4))
  expect_equal(dim(subset(testProtExp, id %in% c('A', 'B'))), c(2, 4))
  expect_equal(dim(subset(testProtExp, id %in% c('A', 'B') & name == 'P1')), c(1, 4))
  expect_equal(dim(subset(testProtExp, id %in% c('A', 'B') & name == 'P3')), c(0, 4))
  expect_error(subset(testProtExp, od %in% c('A', 'B')))

  ## PeptideExperiment
  expect_equal(dim(subset(testPeptExp, id == 'a')), c(1, 4))
  expect_equal(dim(subset(testPeptExp, id %in% c('a', 'b'))), c(2, 4))
  expect_equal(dim(subset(testPeptExp, id %in% c('a', 'b') & name == 'p1')), c(1, 4))
  expect_equal(dim(subset(testPeptExp, id %in% c('a', 'b') & name == 'p3')), c(0, 4))
  expect_error(subset(testPeptExp, od %in% c('A', 'B')))

  ## ProteomicsExperiment
  metaoptions(testPE)[['subsetMode']] <- 'protein'
  metaoptions(testPE)[['linkedSubset']] <- TRUE

  expect_equal(as.vector(dim(subset(testPE, id == 'A'))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B')))), c(2, 3, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B') & name == 'P1'))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B') & name == 'P3'))), c(0, 0, 4, 4))
  expect_error(subset(testPE, od %in% c('A', 'B') & name == 'P3'))

  metaoptions(testPE)[['subsetMode']] <- 'protein'
  metaoptions(testPE)[['linkedSubset']] <- FALSE

  expect_equal(as.vector(dim(subset(testPE, id == 'A'))), c(1, 5, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B')))), c(2, 5, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B') & name == 'P1'))), c(1, 5, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('A', 'B') & name == 'P3'))), c(0, 5, 4, 4))
  expect_error(subset(testPE, od %in% c('A', 'B') & name == 'P3'))

  metaoptions(testPE)[['subsetMode']] <- 'peptide'
  metaoptions(testPE)[['linkedSubset']] <- TRUE

  expect_equal(as.vector(dim(subset(testPE, id == 'a'))), c(1, 1, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b')))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b') & name == 'p1'))), c(1, 1, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b') & name == 'p3'))), c(0, 0, 4, 4))
  expect_error(subset(testPE, od %in% c('a', 'b') & name == 'p3'))

  metaoptions(testPE)[['subsetMode']] <- 'peptide'
  metaoptions(testPE)[['linkedSubset']] <- FALSE

  expect_equal(as.vector(dim(subset(testPE, id == 'a'))), c(3, 1, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b')))), c(3, 2, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b') & name == 'p1'))), c(3, 1, 4, 4))
  expect_equal(as.vector(dim(subset(testPE, id %in% c('a', 'b') & name == 'p3'))), c(3, 0, 4, 4))
  expect_error(subset(testPE, od %in% c('a', 'b') & name == 'p3'))

})

test_that("subsetProt method works", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_equal(dim(subsetProt(testProtExp, id == 'A')), c(1, 4))
  expect_equal(dim(subsetProt(testProtExp, id %in% c('A', 'B'))), c(2, 4))
  expect_equal(dim(subsetProt(testProtExp, id %in% c('A', 'B') & name == 'P1')), c(1, 4))
  expect_equal(dim(subsetProt(testProtExp, id %in% c('A', 'B') & name == 'P3')), c(0, 4))
  expect_error(subsetProt(testProtExp, od %in% c('A', 'B')))

  metaoptions(testPE)[['subsetMode']] <- 'protein'
  metaoptions(testPE)[['linkedSubset']] <- TRUE

  expect_equal(as.vector(dim(subsetProt(testPE, id == 'A'))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B')))), c(2, 3, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B') & name == 'P1'))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B') & name == 'P3'))), c(0, 0, 4, 4))
  expect_error(subsetProt(testPE, od %in% c('A', 'B') & name == 'P3'))

  metaoptions(testPE)[['subsetMode']] <- 'protein'
  metaoptions(testPE)[['linkedSubset']] <- FALSE

  expect_equal(as.vector(dim(subsetProt(testPE, id == 'A'))), c(1, 5, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B')))), c(2, 5, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B') & name == 'P1'))), c(1, 5, 4, 4))
  expect_equal(as.vector(dim(subsetProt(testPE, id %in% c('A', 'B') & name == 'P3'))), c(0, 5, 4, 4))
  expect_error(subsetProt(testPE, od %in% c('A', 'B') & name == 'P3'))

})

test_that("subsetPept method works", {

  ## test objects
  testProtExp <- testList[[1]]
  testPeptExp <- testList[[2]]
  testPE <- testList[[3]]

  expect_equal(dim(subsetPept(testPeptExp, id == 'a')), c(1, 4))
  expect_equal(dim(subsetPept(testPeptExp, id %in% c('a', 'b'))), c(2, 4))
  expect_equal(dim(subsetPept(testPeptExp, id %in% c('a', 'b') & name == 'p1')), c(1, 4))
  expect_equal(dim(subsetPept(testPeptExp, id %in% c('a', 'b') & name == 'p3')), c(0, 4))
  expect_error(subsetPept(testPeptExp, od %in% c('A', 'B')))

  metaoptions(testPE)[['subsetMode']] <- 'peptide'
  metaoptions(testPE)[['linkedSubset']] <- TRUE

  expect_equal(as.vector(dim(subsetPept(testPE, id == 'a'))), c(1, 1, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b')))), c(1, 2, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b') & name == 'p1'))), c(1, 1, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b') & name == 'p3'))), c(0, 0, 4, 4))
  expect_error(subsetPept(testPE, od %in% c('a', 'b') & name == 'p3'))

  metaoptions(testPE)[['subsetMode']] <- 'peptide'
  metaoptions(testPE)[['linkedSubset']] <- FALSE

  expect_equal(as.vector(dim(subsetPept(testPE, id == 'a'))), c(3, 1, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b')))), c(3, 2, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b') & name == 'p1'))), c(3, 1, 4, 4))
  expect_equal(as.vector(dim(subsetPept(testPE, id %in% c('a', 'b') & name == 'p3'))), c(3, 0, 4, 4))
  expect_error(subsetPept(testPE, od %in% c('a', 'b') & name == 'p3'))


})
