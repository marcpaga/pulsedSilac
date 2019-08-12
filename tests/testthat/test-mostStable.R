test_that("mostStable works", {

  PE <- ProtExp(wormsPE)

  expect_error(mostStable(PE, 'asdf'))
  expect_error(stablePE <- mostStable(PE, 'fraction'))
  expect_silent(stablePE <- mostStable(PE, 'fraction', 50))
  expect_equal(nrow(stablePE), 50)

  expect_silent(stablePE <- mostStable(PE, 'fraction', 25))
  expect_equal(nrow(stablePE), 25)

  PE <- ProtExp(wormsPE)[1:20, 1:5]

  assays(PE)[['test']] <- matrix(1:100, nrow = 20, ncol = 5, byrow = TRUE)

  colData(PE)$line <- droplevels(colData(PE)$line)
  expect_silent(stablePE <- mostStable(x = PE,
                                       assayName = 'test',
                                       n = 10))

  expect_equivalent(assays(stablePE)[['test']],
                    assays(PE)[['test']][1:10,])

  PE <- ProtExp(wormsPE)[1:20, c(1:4, 8:11)]

  assays(PE)[['test']] <- cbind(matrix(1:80, nrow = 20, ncol = 4, byrow = TRUE),
                                matrix(c(73:80, 1:72), nrow = 20, ncol = 4,
                                       byrow = TRUE))

  expect_silent(stablePE <- mostStable(x = PE,
                                       assayName = 'test',
                                       n = 3))

  expect_equivalent(assays(stablePE)[['test']],
                    assays(PE)[['test']][3:5,])

})
