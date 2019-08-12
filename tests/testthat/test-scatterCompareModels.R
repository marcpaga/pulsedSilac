test_that("scatterCompareModels works", {

wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')

modelList <- modelTurnover(x = wormsPE[1:10],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t)',
                           start = list(k = 0.02),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = FALSE)

expect_error(scatterCompareModels(modelList = modelList,
                                  conditions = c('OW40', 'OW450'),
                                  value = 'asdf'))

expect_error(scatterCompareModels(modelList = modelList,
                                  conditions = c('OW450'),
                                  value = 'param_values'))

expect_error(scatterCompareModels(modelList = modelList,
                                  conditions = c('asdf','OW450'),
                                  value = 'param_values'))

expect_silent(p <- scatterCompareModels(modelList = modelList,
                                        conditions = c('OW40','OW450'),
                                        value = 'param_values'))

expect_is(p, 'ggplot')
expect_setequal(sapply(p$layers, function(x) class(x$geom)[1]),
                c('GeomPoint', 'GeomAbline'))

expect_silent(p <- scatterCompareModels(modelList = modelList,
                                        conditions = c('OW40','OW450'),
                                        value = 'param_values',
                                        returnDataFrame = TRUE))
expect_is(p, 'data.frame')
expect_named(p, c('OW40', 'OW450', 'param'))
expect_equal(nrow(p), 10)

expect_silent(p <- scatterCompareModels(modelList = modelList,
                                        conditions = c('OW40','OW450'),
                                        value = 'residuals',
                                        returnDataFrame = TRUE))
expect_is(p, 'data.frame')
expect_named(p, c('OW40', 'OW450', 'Time'))
expect_equal(nrow(p), 14)


expect_silent(p <- scatterCompareModels(modelList = modelList,
                                        conditions = c('OW40','OW450'),
                                        value = 'stderror',
                                        returnDataFrame = TRUE))
expect_is(p, 'data.frame')
expect_named(p, c('OW40', 'OW450'))
expect_equal(nrow(p), 10)


modelList <- modelTurnover(x = wormsPE[1:10],
                           assayName = 'fraction',
                           formula = 'fraction ~ 1 - exp(-k*t) + b',
                           start = list(k = 0.02, b = 0),
                           mode = 'protein',
                           robust = FALSE,
                           returnModel = FALSE)

expect_silent(p <- scatterCompareModels(modelList = modelList,
                                        conditions = c('OW40','OW450'),
                                        value = 'param_values',
                                        returnDataFrame = TRUE))
expect_is(p, 'data.frame')
expect_named(p, c('OW40', 'OW450', 'param'))
expect_equal(nrow(p), 20)

})
