test_that("plotDistributionModel works", {

wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')

modelList <- modelTurnover(x = wormsPE[1:30],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

expect_error(plotDistributionModel(modelList = modelList,
                                   value = 'asdf',
                                   plotType = 'density'))

expect_error(plotDistributionModel(modelList = modelList,
                                   value = 'param_values',
                                   plotType = 'asdf'))

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'param_values',
                                         plotType = 'density'))

expect_is(p, 'ggplot')
expect_setequal(sapply(p$layers, function(x) class(x$geom)[1]),
                'GeomDensityRidges')

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'param_values',
                                         plotType = 'density',
                                         returnDataFrame = TRUE))

expect_is(p, 'data.frame')
expect_named(p, c('value', 'condition', 'param'))
expect_equal(nrow(p), 60)
expect_setequal(levels(p$condition), c('OW40', 'OW450'))
expect_setequal(levels(p$param), c('k'))


modelList <- modelTurnover(x = wormsPE[1:30],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t) + b',
                           start = list(k = 0.02, b = 0),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'param_values',
                                         plotType = 'density',
                                         returnDataFrame = TRUE))

expect_is(p, 'data.frame')
expect_named(p, c('value', 'condition', 'param'))
expect_equal(nrow(p), 120)
expect_setequal(levels(p$condition), c('OW40', 'OW450'))
expect_setequal(levels(p$param), c('k', 'b'))


modelList <- modelTurnover(x = wormsPE[1:30, 1:7],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'param_values',
                                         plotType = 'density',
                                         returnDataFrame = TRUE))

expect_is(p, 'data.frame')
expect_named(p, c('value', 'condition', 'param'))
expect_equal(nrow(p), 30)
expect_setequal(levels(p$condition), c('OW40'))
expect_setequal(levels(p$param), c('k'))


modelList <- modelTurnover(x = wormsPE[1:30,],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = TRUE)

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'residuals',
                                         plotType = 'density',
                                         returnDataFrame = TRUE))

expect_is(p, 'data.frame')
expect_named(p, c('value', 'condition', 'time'))
expect_equal(nrow(p), ncol(wormsPE) * 30)
expect_setequal(levels(p$condition), c('OW40', 'OW450'))

expect_silent(p <- plotDistributionModel(modelList = modelList,
                                         value = 'stderror',
                                         plotType = 'density',
                                         returnDataFrame = TRUE))

expect_is(p, 'data.frame')
expect_named(p, c('value', 'condition'))
expect_equal(nrow(p), 60)
expect_setequal(levels(p$condition), c('OW40', 'OW450'))

})
