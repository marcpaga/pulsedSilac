context("test-modelturnover")

test_that("modelturnover proteinExperiment works", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- testPE@ProteinExperiment
  metaoptions(testPE)[['proteinCol']] <- 'Leading.razor.protein'


  expect_silent(ml <- modelTurnover(x = testPE,
                assayName = 'fraction',
                formula = 'fraction ~ 1-exp(-k*t)',
                start = list(k = 0.02),
                robust = FALSE,
                verbose = FALSE,
                returnModel = FALSE))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4)))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), nrow(testPE))
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = TRUE,
                                    verbose = FALSE,
                                    returnModel = FALSE))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror',
                            'weights'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4), 'matrix'))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), nrow(testPE))
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)
  expect_equal(ncol(ml[[7]]), ncol(testPE))

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE))

  expect_is(ml, 'list')
  expect_equal(length(ml), 7)
  expect_equal(length(ml[['models']][[1]]), nrow(testPE))
  expect_equal(length(ml[['models']][[2]]), nrow(testPE))
  expect_is(ml[['models']][[1]][[2]], 'nls')
  expect_equal(names(attributes(ml)), c('names', 'loopCols', 'time', 'cond', 'assayName', 'mode'))

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = TRUE,
                                    verbose = FALSE,
                                    returnModel = TRUE))

  expect_is(ml, 'list')
  expect_equal(length(ml), 8)
  expect_equal(length(ml[['models']][[1]]), nrow(testPE))
  expect_equal(length(ml[['models']][[2]]), nrow(testPE))
  expect_is(ml[['models']][[1]][[2]], 'nls')
  expect_equal(names(attributes(ml)), c('names', 'loopCols', 'time', 'cond', 'assayName', 'mode'))
})


test_that("modelturnover peptideExperiment works", {


  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- testPE@PeptideExperiment
  metaoptions(testPE)[['proteinCol']] <- 'Leading.razor.protein'

  ## 1 model per peptide

  expect_message(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = FALSE,
                                    mode = 'peptide'))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4)))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), nrow(testPE))
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)

  expect_warning(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = TRUE,
                                    verbose = FALSE,
                                    returnModel = FALSE,
                                    mode = 'peptide'))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror',
                            'weights'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4), 'matrix'))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), nrow(testPE))
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)
  expect_equal(ncol(ml[[7]]), ncol(testPE))

  expect_message(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE,
                                    mode = 'peptide'))

  expect_is(ml, 'list')
  expect_equal(length(ml), 7)
  expect_equal(length(ml[['models']][[1]]), nrow(testPE))
  expect_equal(length(ml[['models']][[2]]), nrow(testPE))
  expect_is(ml[['models']][[1]][[6]], 'nls')
  expect_equal(names(attributes(ml)), c('names','loopCols', 'time', 'cond', 'assayName', 'mode'))

  expect_warning(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = TRUE,
                                    verbose = FALSE,
                                    returnModel = TRUE,
                                    mode = 'peptide'))

  expect_is(ml, 'list')
  expect_equal(length(ml), 8)
  expect_equal(length(ml[['models']][[1]]), nrow(testPE))
  expect_equal(length(ml[['models']][[2]]), nrow(testPE))
  expect_is(ml[['models']][[1]][[6]], 'nls')
  expect_equal(names(attributes(ml)), c('names','loopCols', 'time', 'cond', 'assayName', 'mode'))



  ## 1 model per protein

  expect_message(ml <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = FALSE,
                                     mode = 'grouped'))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4)))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), 10)
  expect_equal(nrow(ml[[3]][[1]]), 10)
  expect_equal(nrow(ml[[4]][[1]]), 10)
  expect_equal(nrow(ml[[5]][[1]]), 10)
  expect_equal(nrow(ml[[6]][[1]]), 10)
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)

  expect_warning(ml <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = TRUE,
                                     verbose = FALSE,
                                     returnModel = FALSE,
                                     mode = 'grouped'))

  expect_is(ml, 'list')
  expect_equal(names(ml), c('residuals', 'stderror', 'param_values',
                            'param_pval', 'param_tval', 'param_stderror',
                            'weights'))
  expect_equal(unname(sapply(ml, class)), c(rep('matrix', 2), rep('list', 4), 'matrix'))
  expect_equal(nrow(ml[[1]]), nrow(testPE))
  expect_equal(nrow(ml[[2]]), 10)
  expect_equal(nrow(ml[[3]][[1]]), 10)
  expect_equal(nrow(ml[[4]][[1]]), 10)
  expect_equal(nrow(ml[[5]][[1]]), 10)
  expect_equal(nrow(ml[[6]][[1]]), 10)
  expect_equal(nrow(ml[[7]]), nrow(testPE))
  expect_equal(ncol(ml[[1]]), ncol(testPE))
  expect_equal(ncol(ml[[2]]), 2)
  expect_equal(ncol(ml[[3]][[1]]), 2)
  expect_equal(ncol(ml[[4]][[1]]), 2)
  expect_equal(ncol(ml[[5]][[1]]), 2)
  expect_equal(ncol(ml[[6]][[1]]), 2)
  expect_equal(ncol(ml[[7]]), ncol(testPE))

  expect_message(ml <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = TRUE,
                                     mode = 'grouped'))

  expect_is(ml, 'list')
  expect_equal(length(ml), 7)
  expect_equal(length(ml[['models']][[1]]), 10)
  expect_equal(length(ml[['models']][[2]]), 10)
  expect_is(ml[['models']][[1]][[2]], 'nls')
  expect_equal(names(attributes(ml)), c('names','loopCols', 'time', 'cond', 'prot', 'assayName', 'mode'))


  expect_warning(ml <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = TRUE,
                                     verbose = FALSE,
                                     returnModel = TRUE,
                                     mode = 'grouped'))

  expect_is(ml, 'list')
  expect_equal(length(ml), 8)
  expect_equal(length(ml[['models']][[1]]), 10)
  expect_equal(length(ml[['models']][[2]]), 10)
  expect_is(ml[['models']][[1]][[2]], 'nls')
  expect_equal(names(attributes(ml)), c('names', 'loopCols', 'time', 'cond', 'prot', 'assayName', 'mode'))


})
