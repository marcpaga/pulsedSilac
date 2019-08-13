#' @rdname plotIndividualModel
#' @name plotIndividualModel
#' @title Fitted model(s) for a feature
#'
#' @description Plot the model fit for a specific protein/peptide in different
#' conditions.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object.
#' @param modelList A list containing all the model objects, this should be the
#' output of \code{\link{modelTurnover}} with returnModel as TRUE.
#' @param num The feature number to be plotted.
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#' @param ... Unused.
#'
#' @return A scatter plot with a fitted line or a \code{data.frame}.
#'
#' @examples
#' data('wormsPE')
#' wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
#'
#' modelList <- modelTurnover(x = wormsPE[1:10],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t)',
#'                            start = list(k = 0.02),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' plotIndividualModel(x = wormsPE,
#'                     modelList = modelList,
#'                     num = 2)
#'
#' @import ggplot2
#' @importFrom stats predict
#' @export
setGeneric('plotIndividualModel', function(x, ...){
  standardGeneric('plotIndividualModel')
})

#' @rdname plotIndividualModel
#' @export
setMethod('plotIndividualModel',
          'ProteinExperiment',
          function(x,
                   modelList,
                   num,
                   returnDataFrame = FALSE) {

  ## argument checker ----------------------------------------------------------
  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (!'models' %in% names(modelList)) {
    stop('There are no models in this modelList, did you run modelTurnover ',
         'with returnModel = TRUE?')
  }

  ## data processing -----------------------------------------------------------
  ## models related to that feature
  miniList <- lapply(modelList[['models']], "[[", num)

  ## extract the metaoptions from the attributes of the list
  time <- attributes(modelList)[['time']]
  cond <- attributes(modelList)[['cond']]
  loopCols <- attributes(modelList)[['loopCols']]

  mat <- assaysProt(x)[[attributes(modelList)[['assayName']]]]
  vals <- mat[num, ]

  ## extract the data for each condition and join them in a plotable dataframe
  ## we have two dataframes, one for the original values, and another for the
  ## fitted model. The second one is longer to have a smoother curve.
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      dfOriginalList <- list()
      dfFittedList <- list()
    }
    mod <- miniList[[i]]
    if (!is(mod, 'nls')) {
      warning(paste('No model found for', unique(cond[loopCols[[i]]])))
      next
    }
    origDf <- data.frame(originalval = vals[loopCols[[i]]],
                         time = time[loopCols[[i]]],
                         condition = cond[loopCols[[i]]])
    dfOriginalList[[i]] <- origDf

    all_times <- seq(from = min(time[loopCols[[i]]]),
                     to = max(time[loopCols[[i]]]),
                     by = 1)
    fittedDf <- data.frame(time = all_times,
                           condition = rep(unique(cond[loopCols[[i]]]),
                                           times = length(all_times)),
                           fittedval = predict(mod, list(t = all_times)))
    dfFittedList[[i]] <- fittedDf

  }

  origPlotDf <- do.call('rbind', dfOriginalList)
  fittedPlotDf <- do.call('rbind', dfFittedList)
  origPlotDf$condition <- droplevels(origPlotDf$condition)
  fittedPlotDf$condition <- droplevels(fittedPlotDf$condition)

  if (is.null(origPlotDf)) {
    stop('No models found for any condition')
  }

  ## if the data frames are wanted then they are returned in a list
  if (returnDataFrame) {
    return(list(original_data = origPlotDf, fitted_data = fittedPlotDf))
  }

  ## the plot itself
  ggplot() +
    geom_line(data = fittedPlotDf, aes_string(x = 'time',
                                              y = 'fittedval',
                                              group = 'condition',
                                              color = 'condition')) +
    geom_point(data = origPlotDf, aes_string(x = 'time',
                                             y = 'originalval',
                                             color = 'condition')) +
    scale_color_manual(values = cbPalette) +
    ylab(label = attributes(modelList)[['assayName']]) +
    theme_bw()

})

#' @rdname plotIndividualModel
#' @export
setMethod('plotIndividualModel',
          'PeptideExperiment',
          function(x,
                   modelList,
                   num,
                   returnDataFrame = FALSE) {

  ## argument checker ----------------------------------------------------------
  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (!'models' %in% names(modelList)) {
    stop('There are no models in this modelList, did you run modelTurnover ',
         'with returnModel = TRUE?')
  }

  ## data processing -----------------------------------------------------------
  ## models related to that feature
  miniList <- lapply(modelList[['models']], "[[", num)

  ## extract the metaoptions from the attributes of the list
  time <- attributes(modelList)[['time']]
  cond <- attributes(modelList)[['cond']]
  loopCols <- attributes(modelList)[['loopCols']]

  ## get the rows of the data
  if (attributes(modelList)[['mode']] == 'peptide') {

    mat <- assays(x[num,])[[attributes(modelList)[['assayName']]]]
    repeats <- 1
  } else if (attributes(modelList)[['mode']] == 'grouped') {

    protId <- unique(attributes(modelList)[['prot']])[num]
    rows <- which(attributes(modelList)[['prot']] == protId)
    mat <- assays(x)[[attributes(modelList)[['assayName']]]][rows, , drop = FALSE]
    repeats <- nrow(mat)
  }

  ## extract the data for each condition and join them in a plotable dataframe
  ## we have two dataframes, one for the original values, and another for the
  ## fitted model. The second one is longer to have a smoother curve.
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      dfOriginalList <- list()
      dfFittedList <- list()
    }
    mod <- miniList[[i]]
    if (!is(mod, 'nls')) {
      warning(paste('No model found for', unique(cond[loopCols[[i]]])))
      next
    }

    origDf <- data.frame(originalval = as.vector(mat[,loopCols[[i]]]),
                         time = rep(time[loopCols[[i]]], each = repeats),
                         condition = rep(cond[loopCols[[i]]], each = repeats))

    dfOriginalList[[i]] <- origDf

    all_times <- seq(from = min(time[loopCols[[i]]]),
                     to = max(time[loopCols[[i]]]),
                     by = 1)

    fittedDf <- data.frame(time = all_times,
                           condition = rep(unique(cond[loopCols[[i]]]),
                                           times = length(all_times)),
                           fittedval = predict(mod, list(t = all_times)))

    dfFittedList[[i]] <- fittedDf

  }

  origPlotDf <- do.call('rbind', dfOriginalList)
  fittedPlotDf <- do.call('rbind', dfFittedList)
  origPlotDf$condition <- droplevels(origPlotDf$condition)
  fittedPlotDf$condition <- droplevels(fittedPlotDf$condition)

  if (is.null(origPlotDf)) {
    stop('No models found for any condition')
  }

  ## if the data frames are wanted then they are returned in a list
  if (returnDataFrame) {
    return(list(original_data = origPlotDf, fitted_data = fittedPlotDf))
  }

  ## the plot itself
  ggplot() +
    geom_line(data = fittedPlotDf, aes_string(x = 'time',
                                              y = 'fittedval',
                                              group = 'condition',
                                              color = 'condition')) +
    geom_point(data = origPlotDf, aes_string(x = 'time',
                                             y = 'originalval',
                                             color = 'condition')) +
    scale_color_manual(values = cbPalette) +
    ylab(label = attributes(modelList)[['assayName']]) +
    theme_bw()

})


#' @rdname plotIndividualModel
#' @export
setMethod('plotIndividualModel',
          'ProteomicsExperiment',
          function(x,
                   modelList,
                   num,
                   returnDataFrame = FALSE) {

  if (attributes(modelList)[['mode']] == 'protein') {

    plotIndividualModel(x = x@ProteinExperiment,
                        modelList = modelList,
                        num = num,
                        returnDataFrame = returnDataFrame)

  } else if (attributes(modelList)[['mode']] == 'peptide') {

    plotIndividualModel(x = x@PeptideExperiment,
                        modelList = modelList,
                        num = num,
                        returnDataFrame = returnDataFrame)

  } else if (attributes(modelList)[['mode']] == 'grouped') {

    plotIndividualModel(x = x@PeptideExperiment,
                        modelList = modelList,
                        num = num,
                        returnDataFrame = returnDataFrame)

  }

})
