test_that("plotIndividualModel works", {

wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')

modelList <- modelTurnover(x = wormsPE[1:10],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = FALSE)

expect_error(plotIndividualModel(x = wormsPE,
                                 modelList = modelList,
                                 num = 2))

modelList <- modelTurnover(x = wormsPE[1:10],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

expect_warning(expect_error(plotIndividualModel(x = wormsPE,
                                 modelList = modelList,
                                 num = 1)))

expect_silent(p <- plotIndividualModel(x = wormsPE,
                                       modelList = modelList,
                                       num = 2))

expect_is(p, 'ggplot')
expect_setequal(sapply(p$layers, function(x) class(x$geom)[1]),
             c('GeomLine', 'GeomPoint'))

expect_silent(p <- plotIndividualModel(x = wormsPE,
                                       modelList = modelList,
                                       num = 2,
                                       returnDataFrame = TRUE))

expect_is(p, 'list')
expect_named(p, c('original_data', 'fitted_data'))
expect_named(p[[1]], c('originalval', 'time', 'condition'))
expect_named(p[[2]], c('time', 'condition', 'fittedval'))
expect_setequal(levels(p[[1]]$condition), c('OW40', 'OW450'))
expect_equal(nrow(p[[1]]), ncol(wormsPE))



modelList <- modelTurnover(x = wormsPE[1:10, 1:7],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

p <- plotIndividualModel(x = wormsPE[,1:7],
                         modelList = modelList,
                         num = 2,
                         returnDataFrame = TRUE)

expect_is(p, 'list')
expect_named(p, c('original_data', 'fitted_data'))
expect_named(p[[1]], c('originalval', 'time', 'condition'))
expect_named(p[[2]], c('time', 'condition', 'fittedval'))
expect_setequal(levels(p[[1]]$condition), 'OW40')
expect_equal(nrow(p[[1]]), ncol(wormsPE[,1:7]))

metaoptions(wormsPE)[['proteinCol']] <- 'Protein.group.IDs'
modelList <- modelTurnover(x = wormsPE[1:10,],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'grouped',
                           robust = FALSE,
                           returnModel = TRUE)

p <- plotIndividualModel(x = wormsPE,
                         modelList = modelList,
                         num = 2,
                         returnDataFrame = FALSE)
expect_is(p, 'ggplot')
expect_setequal(sapply(p$layers, function(x) class(x$geom)[1]),
                c('GeomLine', 'GeomPoint'))

p <- plotIndividualModel(x = wormsPE,
                         modelList = modelList,
                         num = 2,
                         returnDataFrame = TRUE)

expect_is(p, 'list')
expect_named(p, c('original_data', 'fitted_data'))
expect_named(p[[1]], c('originalval', 'time', 'condition'))
expect_named(p[[2]], c('time', 'condition', 'fittedval'))
expect_setequal(levels(p[[1]]$condition), c('OW40', 'OW450'))
expect_equivalent(nrow(p[[1]]), nrow(wormsPE[2,])[2] * ncol(wormsPE))



metaoptions(wormsPE)[['proteinCol']] <- 'Protein.group.IDs'
modelList <- modelTurnover(x = wormsPE[1:10,],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'peptide',
                           robust = FALSE,
                           returnModel = TRUE)

p <- plotIndividualModel(x = wormsPE,
                         modelList = modelList,
                         num = 6,
                         returnDataFrame = FALSE)
expect_is(p, 'ggplot')
expect_setequal(sapply(p$layers, function(x) class(x$geom)[1]),
                c('GeomLine', 'GeomPoint'))

p <- plotIndividualModel(x = wormsPE,
                         modelList = modelList,
                         num = 6,
                         returnDataFrame = TRUE)

expect_is(p, 'list')
expect_named(p, c('original_data', 'fitted_data'))
expect_named(p[[1]], c('originalval', 'time', 'condition'))
expect_named(p[[2]], c('time', 'condition', 'fittedval'))
expect_setequal(levels(p[[1]]$condition), c('OW40', 'OW450'))
expect_equivalent(nrow(p[[1]]), ncol(wormsPE))


})
