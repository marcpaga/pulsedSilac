#' @rdname plotDistributionAssay
#' @name plotDistributionAssay
#' @title Distribution of assay data per condition and timepoint.
#'
#' @description Plot the distribution of the data stored in an assay using
#' boxplots or density distributions.
#'
#' @param x A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or a
#' \code{SilacProteomicsExperiment} object.
#' @param assayName Name of the assay to use in the plot.
#' @param plotType  A \code{character} indicating which geometry to plot:
#' 'boxplot' or 'density'. (default = 'density')
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#' @param mode A \code{character} indicating which level of data to use,
#' either "protein" or "peptide". Only relevant for ProteomicsExperiment
#' inputs.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param timeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different timepoints.
#' @param ... Unused.
#'
#' @return A ggplot2 object or a \code{data.frame} with the data that would be
#' plotted.
#'
#' @examples
#' data('wormsPE')
#' plotDistributionAssay(wormsPE, assayName = 'ratio')
#'
#' @importFrom ggridges geom_density_ridges
#' @import ggplot2
#' @export
setGeneric('plotDistributionAssay', function(x, ...){
  standardGeneric('plotDistributionAssay')
})

#' @rdname plotDistributionAssay
#' @export
setMethod('plotDistributionAssay', 'SilacProteinExperiment',
          function(x,
                   assayName,
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  ## argument checker ----------------------------------------------------------
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  if (!plotType %in% c('boxplot', 'density')) {
    txt <- c('plotType must be "boxplot" or "density"')
    stop(txt)
  }

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (!missing(conditionCol)) {
    metadata(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metadata(x)[['timeCol']] <- timeCol
  }

  ## Data and options processing -----------------------------------------------
  mat <- assays(x)[[assayName]]

  loopCols <- .loopWrapper(x, 'conditionCol')

  ## make a data.frame for each condition
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      dfList <- list()
    }

    timeCol <- .giveMetaoption(x, 'timeCol')
    if (is.na(timeCol)) {
      timeVec <- seq_along(loopCols[[i]])
    } else {
      timeVec <- colData(x)[, .giveMetaoption(x, 'timeCol')][loopCols[[i]]]
    }

    values <- as.vector(mat[, loopCols[[i]]])
    tempDf <- data.frame(value = values,
                         time = rep(timeVec, each = nrow(mat)))
    if (!is.null(names(loopCols)[i])) {
      tempDf$condition <- names(loopCols)[i]
    } else {
      tempDf$condition <- 'condition'
    }


    dfList[[i]] <- tempDf
  }

  ## join all the data.frames
  plotDf <- do.call('rbind', dfList)
  plotDf$time <- as.factor(plotDf$time)
  plotDf$condition <- as.factor(plotDf$condition)


  if (returnDataFrame) {
    colnames(plotDf)[1] <- assayName
    return(plotDf)
  }

  ## plotting ------------------------------------------------------------------
  if (plotType == 'density') {

    p <- ggplot(data = plotDf) +
      geom_density_ridges(aes_string(x = 'value',
                                     y = 'time',
                                     fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~condition) +
      labs(x = assayName) +
      theme_bw()

  } else if (plotType == 'boxplot') {

    p <- ggplot(data = plotDf) +
      geom_boxplot(aes_string(x = 'time',
                              y = 'value',
                              fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      labs(y = assayName) +
      theme_bw()

  }

  p

})

#' @rdname plotDistributionAssay
#' @export
setMethod('plotDistributionAssay', 'SilacPeptideExperiment',
          function(x,
                   assayName,
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  callNextMethod()

})

#' @rdname plotDistributionAssay
#' @export
setMethod('plotDistributionAssay', 'SilacProteomicsExperiment',
          function(x,
                   assayName,
                   mode = 'protein',
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  experiment <- switch(mode,
                       protein = x@SilacProteinExperiment,
                       peptide = x@SilacPeptideExperiment)

  plotDistributionAssay(x = experiment,
            assayName = assayName,
            plotType = plotType,
            returnDataFrame = returnDataFrame,
            conditionCol = conditionCol,
            timeCol = timeCol)


})
