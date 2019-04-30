context("test-proteomicsexperiment")

test_that("constructor works", {

  assay1 <- matrix(1:9, ncol = 3, nrow = 3)
  assay2 <- matrix(9:1, ncol = 3, nrow = 3)
  assays_list <- list(a1 = assay1, a2 = assay2)
  colData <- data.frame(sample = LETTERS[1:3])
  rowData <- data.frame(peptide = letters[1:3])

  proExp <- ProteinExperiment(assays = assays_list,
                              rowData = rowData,
                              colData = colData)

  assay1 <- matrix(1:12, ncol = 3, nrow = 4)
  assay2 <- matrix(12:1, ncol = 3, nrow = 4)
  assays_list <- list(a1 = assay1, a2 = assay2)
  colData <- data.frame(sample = LETTERS[1:3])
  rowData <- data.frame(peptide = letters[1:4])

  pepExp <- PeptideExperiment(assays = assays_list,
                              rowData = rowData,
                              colData = colData)

  expect_silent(PE <- ProteomicsExperiment(ProteinExperiment = proExp,
                                           PeptideExperiment = pepExp))
  expect_equal(length(metadata(PE)), 0)
  expect_equal(length(metaoptions(PE)), 9)


})
