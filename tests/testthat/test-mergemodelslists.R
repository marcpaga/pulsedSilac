context("test-mergemodelslists")

test_that("mergemodelslists works protein", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  metaoptions(testPE)[['proteinCol']] <- 'Leading.razor.protein'

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'protein',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))
})


test_that("mergemodelslists works grouped", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  metaoptions(testPE)[['proteinCol']] <- 'Leading.razor.protein'

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'grouped',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))
})



test_that("mergemodelslists works peptide", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  metaoptions(testPE)[['proteinCol']] <- 'Leading.razor.protein'

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = FALSE)

  ml3 <- mergeModelsLists(ml1, ml2)
  expect_identical(ml0, ml3)

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = TRUE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))

  ml0 <- modelTurnover(x = testPE,
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml1 <- modelTurnover(x = testPE[,1:7],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml2 <- modelTurnover(x = testPE[,8:14],
                       assayName = 'fraction',
                       formula = 'fraction ~ 1-exp(-k*t)',
                       mode = 'peptide',
                       start = list(k = 0.02),
                       robust = FALSE,
                       returnModel = TRUE)

  ml3 <- mergeModelsLists(ml1, ml2)
  ## cannot use identical here because the model objects have some differences
  expect_equal(ml0, ml3)
  expect_identical(attributes(ml0), attributes(ml3))
})
