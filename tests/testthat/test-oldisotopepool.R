test_that("addMiscleavedPeptides works", {

  PE <- testList[[1]]
  recycleDf <- data.frame(id = rep(letters[1:4], each = 2),
                          mod = rep(c('NN', 'ON'), times = 4),
                          int1 = c(1:8),
                          int2 = c(1:8),
                          int3 = c(1:8),
                          int4 = c(1:8),
                          oth1 = letters[16:23],
                          oth2 = letters[18:25])

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'mod',
                                                dataCols = 3:7))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'asdf',
                                                dataCols = 3:6))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'asdf',
                                                modCol = 'mod',
                                                dataCols = 3:6))

  expect_silent(missPE <- addMisscleavedPeptides(x = PE,
                                                 newdata = recycleDf,
                                                 idColPept = 'id',
                                                 modCol = 'mod',
                                                 dataCols = 3:6))

  expect_is(missPE, 'PeptideExperiment')
  expect_equal(dim(missPE), c(4, 4))
  expect_length(assays(missPE), 2)
  expect_named(assays(missPE), c('NN', 'ON'))



  PE <- testList[[2]]
  recycleDf <- data.frame(id = rep(letters[1:4], each = 2),
                          mod = rep(c('NN', 'ON'), times = 4),
                          int1 = c(1:8),
                          int2 = c(1:8),
                          int3 = c(1:8),
                          int4 = c(1:8),
                          oth1 = letters[16:23],
                          oth2 = letters[18:25])

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'mod',
                                                dataCols = 3:7))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'asdf',
                                                dataCols = 3:6))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'asdf',
                                                modCol = 'mod',
                                                dataCols = 3:6))

  expect_silent(missPE <- addMisscleavedPeptides(x = PE,
                                                 newdata = recycleDf,
                                                 idColPept = 'id',
                                                 modCol = 'mod',
                                                 dataCols = 3:6))

  expect_is(missPE, 'PeptideExperiment')
  expect_equal(dim(missPE), c(5, 4))
  expect_length(assays(missPE), 4)
  expect_named(assays(missPE)[3:4], c('NN', 'ON'))


  PE <- testList[[3]]
  recycleDf <- data.frame(id = rep(letters[1:4], each = 2),
                          mod = rep(c('NN', 'ON'), times = 4),
                          int1 = c(1:8),
                          int2 = c(1:8),
                          int3 = c(1:8),
                          int4 = c(1:8),
                          oth1 = letters[16:23],
                          oth2 = letters[18:25])

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'mod',
                                                dataCols = 3:7))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'id',
                                                modCol = 'asdf',
                                                dataCols = 3:6))

  expect_error(missPE <- addMisscleavedPeptides(x = PE,
                                                newdata = recycleDf,
                                                idColPept = 'asdf',
                                                modCol = 'mod',
                                                dataCols = 3:6))

  expect_silent(missPE <- addMisscleavedPeptides(x = PE,
                                                 newdata = recycleDf,
                                                 idColPept = 'id',
                                                 modCol = 'mod',
                                                 dataCols = 3:6))

  expect_is(missPE, 'ProteomicsExperiment')
  expect_equal(unname(dim(missPE)[2, ]), c(5, 4))
  expect_length(assaysPept(missPE), 4)
  expect_named(assaysPept(missPE)[3:4], c('NN', 'ON'))

})

test_that("calculateOldIsotopePool works", {



})
