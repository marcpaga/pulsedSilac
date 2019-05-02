context("test-rbind")

test_that("rbind works", {

  ## ProteinExperiment
  testPE <- testList[[1]]

  expect_silent(testPE2 <- rbind(testPE[1:2,], testPE[3,]))
  expect_equivalent(testPE, testPE2)
  expect_equivalent(colData(testPE), colData(testPE2))
  expect_equivalent(rowData(testPE), rowData(testPE2))
  expect_equivalent(metadata(testPE), metadata(testPE2))
  expect_equivalent(metaoptions(testPE), metaoptions(testPE2))

  ## PeptideExperiment
  testPE <- testList[[2]]

  expect_silent(testPE2 <- rbind(testPE[1:2,], testPE[3:5,]))
  expect_equivalent(testPE, testPE2)
  expect_equivalent(colData(testPE), colData(testPE2))
  expect_equivalent(rowData(testPE), rowData(testPE2))
  expect_equivalent(metadata(testPE), metadata(testPE2))
  expect_equivalent(metaoptions(testPE), metaoptions(testPE2))

  ## ProteomicsExperiment
  testPE <- testList[[3]]

  expect_silent(testPE2 <- rbind(testPE[1:2,], testPE[3,]))
  expect_equivalent(testPE@ProteinExperiment, testPE2@ProteinExperiment)
  expect_equivalent(testPE@PeptideExperiment, testPE2@PeptideExperiment)
  expect_equivalent(colData(testPE), colData(testPE2))
  expect_equivalent(rowDataProt(testPE), rowDataProt(testPE2))
  expect_equivalent(rowDataPept(testPE), rowDataPept(testPE2))
  expect_equivalent(metadata(testPE), metadata(testPE2))
  expect_equivalent(metaoptions(testPE), metaoptions(testPE2))
  expect_equivalent(linkerDf(testPE), linkerDf(testPE2))



})
