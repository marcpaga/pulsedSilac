context("test-proteinexperiment")

test_that("constructor works", {

  ## the basic
  assay1 <- matrix(1:9, ncol = 3, nrow = 3)
  assay2 <- matrix(9:1, ncol = 3, nrow = 3)
  assays_list <- list(a1 = assay1, a2 = assay2)
  colData <- data.frame(sample = LETTERS[1:3])
  rowData <- data.frame(peptide = letters[1:3])

  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData))

  expect_equal(length(assays(proExp)), 2)
  expect_equal(nrow(rowData(proExp)), 3)
  expect_equal(ncol(rowData(proExp)), 1)
  expect_equal(nrow(colData(proExp)), 3)
  expect_equal(ncol(colData(proExp)), 1)
  expect_equal(length(metadata(proExp)), 0)
  expect_equal(length(metaoptions(proExp)), 2)

  ## without assays
  expect_error(proExp <- ProteinExperiment(rowData = rowData,
                                           colData = colData))

  ## without rowData
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            colData = colData))

  ## without colData
  expect_error(proExp <- ProteinExperiment(assays = assays_list,
                                           rowData = rowData))

  ## with added metadata
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData,
                                            metadata = list(author = 'me')))
  expect_equal(length(metadata(proExp)), 1)
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData,
                                            metadata = list(author = 'me'),
                                            conditionCol = 'sample'))
  expect_equal(length(metadata(proExp)), 1)
  expect_equal(length(metaoptions(proExp)), 2)

  colData <- data.frame(sample = LETTERS[1:3],
                        time = c(1:3),
                        repInt = c(1:3),
                        repTime = c(1, 1, 1))

  ## all metadata given
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData,
                                            metadata = list(author = 'me'),
                                            conditionCol = 'sample',
                                            timeCol = 'time'))

  ## mixed character and numeric metadata
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData,
                                            metadata = list(author = 'me'),
                                            conditionCol = 'sample',
                                            timeCol = 2))

  ## metaoptions target the same colData columns
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                           rowData = rowData,
                                           colData = colData,
                                           metadata = list(author = 'me'),
                                           conditionCol = 2,
                                           timeCol = 2))


  ## metaoptions in both direct argument and metadata argument
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                           rowData = rowData,
                                           colData = colData,
                                           metadata = list(author = 'me', timeCol = 2),
                                           conditionCol = 1))

  ## user metaoptions
  expect_silent(proExp <- ProteinExperiment(assays = assays_list,
                                            rowData = rowData,
                                            colData = colData,
                                            metadata = list(author = 'me', timeCol = 2),
                                            conditionCol = 1,
                                            metaoptions = list(a = 1)))
  expect_length(metaoptions(proExp), 3)
  expect_named(metaoptions(proExp), c('conditionCol', 'timeCol', 'a'), ignore.order = TRUE)

  expect_error(proExp <- ProteinExperiment(assays = assays_list,
                                           rowData = rowData,
                                           colData = colData,
                                           metadata = list(author = 'me', timeCol = 2),
                                           conditionCol = 1,
                                           metaoptions = list(timeCol = 1)))

})

