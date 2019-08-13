test_that("calculateAIC works", {
  data('wormsPE')
  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- ProtExp(testPE)


  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = FALSE))

  expect_error(calculateAIC(data.frame()))
  expect_error(calculateAIC(ml))

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE))

  expect_silent(ml <- calculateAIC(ml))
  expect_length(ml, 8)
  expect_named(ml[8], 'AIC')
  expect_is(ml[[8]], 'matrix')

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE))

  expect_silent(ml <- calculateAIC(ml, smallSampleSize = TRUE))
  expect_length(ml, 8)
  expect_named(ml[8], 'AIC')
  expect_is(ml[[8]], 'matrix')


})


test_that("compareAIC works", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- ProtExp(testPE)


  expect_silent(ml1 <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = TRUE))
  expect_silent(ml1 <- calculateAIC(ml1))

  expect_error(compareAIC(ml1))

  expect_silent(ml2 <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t) + b',
                                     start = list(k = 0.02, b = 0),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = TRUE))

  expect_error(compareAIC(ml1, ml2))
  expect_silent(ml2 <- calculateAIC(ml2))

  expect_silent(ml_compare <- compareAIC(ml1, ml2))
  expect_is(ml_compare, 'list')
  expect_length(ml_compare, 1)
  expect_length(unlist(ml_compare, FALSE), 2)

  expect_silent(ml3 <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1.2-exp(-k*t) + b',
                                     start = list(k = 0.02, b = 0),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = TRUE))
  expect_silent(ml3 <- calculateAIC(ml3))

  expect_silent(ml_compare <- compareAIC(ml1, ml2, ml3))
  expect_is(ml_compare, 'list')
  expect_length(ml_compare, 1)
  expect_length(unlist(ml_compare, FALSE), 3)




})
