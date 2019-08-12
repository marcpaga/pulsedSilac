test_that("calculateIsotopeRatio works", {

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

  expect_error(calculateIsotopeRatio(PE, newIsotopeAssay = 'new',
                                     oldIsotopeAssay = 'asdf'))

  expect_error(calculateIsotopeRatio(PE, newIsotopeAssay = 'asdf',
                                     oldIsotopeAssay = 'old'))

  expect_silent(PE <- calculateIsotopeRatio(PE, newIsotopeAssay = 'new',
                                      oldIsotopeAssay = 'old'))

  expect_length(assays(PE), 4)
  expect_named(assays(PE)[2], 'ratio')
  expect_equivalent(assays(PE)[['ratio']],
                    matrix(c(NA, 1, 1, 2,
                             1/2, 1, 1, 2,
                             0, 1, 1, NA),
                           nrow = 3, ncol = 4, byrow = TRUE))



})
