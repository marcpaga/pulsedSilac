test_that("calculateAIC works", {

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

})
