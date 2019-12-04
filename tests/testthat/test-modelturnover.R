context("test-modelturnover")

test_that("modelturnover proteinExperiment works", {

  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- testPE@SilacProteinExperiment
  metadata(testPE)[['proteinCol']] <- 'Leading.razor.protein'


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

  colnames(testPE) <- LETTERS[1:14]
  rownames(testPE) <- LETTERS[1:10]

  expect_silent(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE))
  expect_equal(rownames(ml[[1]]), rownames(testPE))
  expect_equal(colnames(ml[[1]]), colnames(testPE))
  expect_equal(rownames(ml[[2]]), rownames(testPE))
  expect_equal(colnames(ml[[2]]), c('OW40', 'OW450'))

})


test_that("modelturnover peptideExperiment works", {


  wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
  testPE <- wormsPE[1:10,]
  testPE <- testPE@SilacPeptideExperiment
  metadata(testPE)[['proteinCol']] <- 'Leading.razor.protein'

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

  colnames(testPE) <- LETTERS[1:14]
  rownames(testPE) <- rowData(testPE)$Sequence

  expect_message(ml <- modelTurnover(x = testPE,
                                    assayName = 'fraction',
                                    formula = 'fraction ~ 1-exp(-k*t)',
                                    start = list(k = 0.02),
                                    robust = FALSE,
                                    verbose = FALSE,
                                    returnModel = TRUE,
                                    mode = 'peptide'))
  expect_equal(rownames(ml[[1]]), rownames(testPE))
  expect_equal(colnames(ml[[1]]), colnames(testPE))
  expect_equal(rownames(ml[[2]]), rownames(testPE))
  expect_equal(colnames(ml[[2]]), c('OW40', 'OW450'))


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


  expect_message(ml <- modelTurnover(x = testPE,
                                     assayName = 'fraction',
                                     formula = 'fraction ~ 1-exp(-k*t)',
                                     start = list(k = 0.02),
                                     robust = FALSE,
                                     verbose = FALSE,
                                     returnModel = TRUE,
                                     mode = 'grouped'))
  expect_equal(rownames(ml[[1]]), rownames(testPE))
  expect_equal(colnames(ml[[1]]), colnames(testPE))
  expect_equal(rownames(ml[[2]]),
               unique(rowData(testPE)[,metaoptions(testPE)[['proteinCol']]]))
  expect_equal(colnames(ml[[2]]), c('OW40', 'OW450'))


})


test_that("modelturnover differet nls algorithms works", {


  frac <- 1 - exp(-0.02 * c(4, 8, 16, 24, 36, 48)) + rnorm(6, mean = 0, sd = 0.02)
  data <- data.frame(t = c(4, 8, 16, 24, 36, 48),
                     fraction =frac)
  formula <- as.formula('fraction ~ 1 - exp(-k * t)')
  start <- list(k = 0.025)
  robust <- FALSE

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE)

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE,
                              algorithm = 'plinear')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror, length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))


  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE,
                              algorithm = 'port')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))


  frac <- 1 - exp(-0.02 * c(4, 8, 16, 24, 36, 48)) + rnorm(6, mean = 0, sd = 0.02)
  data <- data.frame(t = c(4, 8, 16, 24, 36, 48),
                     fraction =frac)
  formula <- as.formula('fraction ~ 1 - exp(-k * t) + b')
  start <- list(k = 0.025, b = 0)
  robust <- FALSE

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE)

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE,
                              algorithm = 'plinear')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror, length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))


  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE,
                              algorithm = 'port')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))



})



test_that("modelturnover differet nlrob algorithms works", {


  frac <- 1 - exp(-0.02 * c(4, 8, 16, 24, 36, 48)) + rnorm(6, mean = 0, sd = 0.02)
  data <- data.frame(t = c(4, 8, 16, 24, 36, 48),
                     fraction =frac)
  formula <- as.formula('fraction ~ 1 - exp(-k * t) + b')
  start <- list(k = 0.025, b = 0)
  lower <- c(k = 0, b = 0)
  upper <- c(k = 0.5, b = 0.1)
  robust <- TRUE

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              start = start,
                              robust = robust,
                              returnModel = FALSE,
                              method = 'M')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))

  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              lower = lower,
                              upper = upper,
                              robust = robust,
                              returnModel = FALSE,
                              method = 'MM')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror, length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))


  modeldata <- .modelTurnover(data = data,
                              formula = formula,
                              lower = lower,
                              upper = upper,
                              robust = robust,
                              returnModel = FALSE,
                              method = 'tau')

  expect_length(modeldata$residuals, nrow(data))
  expect_length(modeldata$stderror, 1)
  expect_length(modeldata$params.vals, length(start))
  expect_length(modeldata$params.stderror,length(start))
  expect_length(modeldata$params.tval, length(start))
  expect_length(modeldata$params.pval, length(start))

  modeldata <- expect_error(.modelTurnover(data = data,
                              formula = formula,
                              lower = lower,
                              upper = upper,
                              robust = robust,
                              returnModel = FALSE,
                              method = 'CM'))


  modeldata <- expect_error(.modelTurnover(data = data,
                              formula = formula,
                              lower = lower,
                              upper = upper,
                              robust = robust,
                              returnModel = FALSE,
                              method = 'mtl', silent = FALSE))
})


