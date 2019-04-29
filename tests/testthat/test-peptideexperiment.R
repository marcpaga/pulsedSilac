context("test-peptideexperiment")

test_that("constructor works", {

## the basic
assay1 <- matrix(1:12, ncol = 3, nrow = 4)
assay2 <- matrix(12:1, ncol = 3, nrow = 4)
assays_list <- list(a1 = assay1, a2 = assay2)
colData <- data.frame(sample = LETTERS[1:3])
rowData <- data.frame(peptide = letters[1:4])

expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData))

expect_equal(length(assays(pepExp)), 2)
expect_equal(nrow(rowData(pepExp)), 4)
expect_equal(ncol(rowData(pepExp)), 1)
expect_equal(nrow(colData(pepExp)), 3)
expect_equal(ncol(colData(pepExp)), 1)
expect_equal(length(metadata(pepExp)), 4)

## without assays
expect_error(pepExp <- PeptideExperiment(rowData = rowData,
                                         colData = colData))

## without rowData
expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                           colData = colData))

## without colData
expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData))

## with added metadata
expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData,
                                          metadata = list(author = 'me')))
expect_equal(length(metadata(pepExp)), 5)
expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData,
                                          metadata = list(author = 'me'),
                                          conditionCol = 'sample'))
expect_equal(length(metadata(pepExp)), 5)
expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData,
                                          metadata = list(author = 'me'),
                                          conditionCol = '1',
                                          timeCol = '2',
                                          replicateIntCol = '3',
                                          replicateTimeCol = '4'))

expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = list(author = 'me'),
                                         replicateTimeCol = 5))

colData <- data.frame(sample = LETTERS[1:3],
                      time = c(1:3),
                      repInt = c(1:3),
                      repTime = c(1, 1, 1))

## all metadata given
expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = list(author = 'me'),
                                         conditionCol = 'sample',
                                         timeCol = 'time',
                                         replicateIntCol = 'repInt',
                                         replicateTimeCol = 'repTime'))

## mixed character and numeric metadata
expect_silent(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData,
                                          metadata = list(author = 'me'),
                                          conditionCol = 'sample',
                                          timeCol = 2,
                                          replicateIntCol = 'repInt',
                                          replicateTimeCol = 4))

## metaoptions target the same colData columns
expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                          rowData = rowData,
                                          colData = colData,
                                          metadata = list(author = 'me'),
                                          conditionCol = 2,
                                          timeCol = 2,
                                          replicateIntCol = 'repInt',
                                          replicateTimeCol = 4))

expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = list(author = 'me'),
                                         conditionCol = 'repInt',
                                         timeCol = 2,
                                         replicateIntCol = 'repInt',
                                         replicateTimeCol = 4))

expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = list(author = 'me'),
                                         conditionCol = 'repInt',
                                         timeCol = 'repInt',
                                         replicateIntCol = 'repInt',
                                         replicateTimeCol = 'repInt'))

expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                                         rowData = rowData,
                                         colData = colData,
                                         metadata = list(author = 'me'),
                                         conditionCol = 1,
                                         timeCol = 'repInt',
                                         replicateIntCol = 1,
                                         replicateTimeCol = 'repInt'))

## metaoptions in both direct argument and metadata argument

expect_error(pepExp <- PeptideExperiment(assays = assays_list,
                            rowData = rowData,
                            colData = colData,
                            metadata = list(author = 'me', timeCol = 2),
                            conditionCol = 1))
})
