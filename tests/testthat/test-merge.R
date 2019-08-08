context("test-merge")

test_that("merge ProteinExperiment", {

  testPE <- testList[[1]]

  testPE2 <- testPE
  rowData(testPE2) <- data.frame(id = LETTERS[2:4],
                                 name = c('P2', 'P3', 'P4'))
  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equal(nrow(testPE3), 4)
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assays(testPE3)), length(unique(c(assayNames(testPE),
                                                        assayNames(testPE2)))))

  testPE <- testList[[1]]

  testPE2 <- testPE
  rowData(testPE2) <- data.frame(id = LETTERS[2:4],
                                 name = c('P2', 'P3', 'P4'))
  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  names(assays(testPE)) <- NULL
  names(assays(testPE2)) <- NULL
  assays(testPE2)[[3]] <- assays(testPE2)[[2]]

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equal(nrow(testPE3), 4)
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assays(testPE3)), max(length(assays(testPE)),
                                            length(assays(testPE2))))


})

test_that("merge PeptideExperiment", {

  testPE <- testList[[2]]

  testPE2 <- testPE
  rowData(testPE2) <- data.frame(id = letters[2:6],
                                 name = c('p2', 'p3', 'p4', 'p5', 'p6'))
  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equal(nrow(testPE3), 6)
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assays(testPE3)), length(unique(c(assayNames(testPE),
                                                        assayNames(testPE2)))))

  testPE <- testList[[2]]

  testPE2 <- testPE
  rowData(testPE2) <- data.frame(id = letters[2:6],
                                 name = c('p2', 'p3', 'p4', 'p5', 'p6'))
  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  names(assays(testPE)) <- NULL
  names(assays(testPE2)) <- NULL
  assays(testPE2)[[3]] <- assays(testPE2)[[2]]

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equal(nrow(testPE3), 6)
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assays(testPE3)), max(length(assays(testPE)),
                                            length(assays(testPE2))))


})

test_that("merge ProteomicsExperiment", {

  testPE <- testList[[3]]

  testPE2 <- testPE
  rowDataProt(testPE2) <- data.frame(id = LETTERS[2:4],
                                     name = c('P2', 'P3', 'P4'))
  rowDataPept(testPE2) <- data.frame(id = letters[2:6],
                                     name = c('p2', 'p3', 'p4', 'p5', 'p6'))

  linkerDf <- buildLinkerDf(protIDs = LETTERS[2:4],
                            pepIDs = letters[2:6],
                            protToPep = list(B = c('b','c'),
                                             C = c('d', 'e'),
                                             D = c('f')))
  linkerDf(testPE2) <- linkerDf

  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equivalent(nrow(testPE3), c(4, 6))
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assaysProt(testPE3)), length(unique(c(assayNamesProt(testPE),
                                                        assayNamesProt(testPE2)))))
  expect_equal(length(assaysPept(testPE3)), length(unique(c(assayNamesPept(testPE),
                                                            assayNamesPept(testPE2)))))

  testPE <- testList[[3]]

  testPE2 <- testPE
  rowDataProt(testPE2) <- data.frame(id = LETTERS[2:4],
                                     name = c('P2', 'P3', 'P4'))

  rowDataPept(testPE2) <- data.frame(id = letters[2:6],
                                     name = c('p2', 'p3', 'p4', 'p5', 'p6'))
  colData(testPE2) <- data.frame(sample = 1:4,
                                 condition = c(1, 1, 2, 2),
                                 time = c(3, 4, 3, 4))

  linkerDf <- buildLinkerDf(protIDs = LETTERS[2:4],
                            pepIDs = letters[2:6],
                            protToPep = list(B = c('b','c'),
                                             C = c('d', 'e'),
                                             D = c('f')))
  linkerDf(testPE2) <- linkerDf

  names(assaysProt(testPE)) <- NULL
  names(assaysPept(testPE)) <- NULL
  names(assaysProt(testPE2)) <- NULL
  names(assaysPept(testPE2)) <- NULL
  assaysProt(testPE2)[[3]] <- assaysProt(testPE2)[[2]]
  assaysPept(testPE2)[[3]] <- assaysPept(testPE2)[[2]]

  testPE3 <- merge(testPE, testPE2)
  expect_equal(ncol(testPE3), ncol(testPE) + ncol(testPE2))
  expect_equivalent(nrow(testPE3), c(4, 6))
  expect_equal(metadata(testPE3), metadata(testPE))
  expect_equal(metaoptions(testPE3), metaoptions(testPE))
  expect_equal(length(assays(testPE3)), max(length(assays(testPE)),
                                            length(assays(testPE2))))


})
