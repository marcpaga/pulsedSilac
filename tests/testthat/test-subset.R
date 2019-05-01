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

  protExp <- ProteinExperiment(assays = protmatrix,
                               rowData = protdata,
                               colData = coldata)
  peptExp <- PeptideExperiment(assays = pepmatrix,
                               rowData = pepdata,
                               colData = coldata)


  PE <- ProteomicsExperiment(ProteinExperiment = protExp,
                             PeptideExperiment = peptExp,
                             colData = coldata,
                             linkerDf = linkM,
                             idColProt = 'a',
                             idColPept = 'b')

  metadata(PE) <- list(a = 1, b = 2)

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
