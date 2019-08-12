test_that("calculateIsotopeFraction simple", {

  PE <- testList[[1]]
  expect_error(PE <- calculateIsotopeFraction(PE, 'asdf'))
  expect_silent(PE <- calculateIsotopeFraction(PE, 'ratio'))
  expect_length(assays(PE), 3)
  expect_named(assays(PE)[3], 'fraction')

  expect_equal(assays(PE)[['fraction']],
               assays(PE)[['ratio']]/(assays(PE)[['ratio']] + 1))


})

test_that("calculateIsotopeFraction complex", {

  PE <- testList[[1]]
  assays(PE)[['new']] <- matrix(data = c(NA, 1, 1, 2,
                                         1, 1, 1, 2,
                                         0, 1, 1, 2),
                                nrow = 3,
                                ncol = 4, byrow = TRUE)
  assays(PE)[['old']] <- matrix(data = c(2, 1, 1, 1,
                                         2, 1, 1, 1,
                                         2, 1, 1, NA),
                                nrow = 3,
                                ncol = 4, byrow = TRUE)

  expect_error(calculateIsotopeFraction(PE,
                                        oldIsoAssay = 'old',
                                        newIsoAssay = 'new',
                                        earlyTimepoints,
                                        lateTimepoints))

  expect_error(calculateIsotopeFraction(PE,
                                        oldIsoAssay = 'old1',
                                        newIsoAssay = 'new',
                                        earlyTimepoints = 1,
                                        lateTimepoints))

  expect_error(calculateIsotopeFraction(PE,
                                        oldIsoAssay = 'old',
                                        newIsoAssay = 'new1',
                                        earlyTimepoints = 1,
                                        lateTimepoints))

  expect_silent(PE <- calculateIsotopeFraction(PE,
                                               oldIsoAssay = 'old',
                                               newIsoAssay = 'new',
                                               earlyTimepoints = 1,
                                               lateTimepoints))

  expect_length(assays(PE), 5)
  expect_named(assays(PE)[5], 'fraction')

  expect_equivalent(assays(PE)[['fraction']],
                    matrix(c(0, 0.5, 0.5, 2/3,
                             1/3, 0.5, 0.5, 2/3,
                             0, 0.5, 0.5, NA),
                      nrow = 3, ncol = 4, byrow = TRUE))

  expect_silent(PE <- calculateIsotopeFraction(x = PE,
                                               oldIsoAssay = 'old',
                                               newIsoAssay = 'new',
                                               lateTimepoints = 4))

  expect_equivalent(assays(PE)[['fraction']],
                    matrix(c(NA, 0.5, 0.5, 2/3,
                             1/3, 0.5, 0.5, 2/3,
                             0, 0.5, 0.5, 1),
                           nrow = 3, ncol = 4, byrow = TRUE))

  expect_silent(PE <- calculateIsotopeFraction(x = PE,
                                               oldIsoAssay = 'old',
                                               newIsoAssay = 'new',
                                               earlyTimepoints = 1,
                                               lateTimepoints = 4))

  expect_equivalent(assays(PE)[['fraction']],
                    matrix(c(0, 0.5, 0.5, 2/3,
                             1/3, 0.5, 0.5, 2/3,
                             0, 0.5, 0.5, 1),
                           nrow = 3, ncol = 4, byrow = TRUE))

})
