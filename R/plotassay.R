#' @rdname boxplotAssay
#' @name boxplotAssay
#' @title Distribution of assay data per condition and timepoint.
#'
#' @description Plot the distribution of the data stored in an assay using
#' boxplots or density distributions.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
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
#' boxplotAssay(wormsPE, assayName = 'ratio')
#'
#' @importFrom ggridges geom_density_ridges
#' @import ggplot2
#' @export
setGeneric('boxplotAssay', function(x, ...){
  standardGeneric('boxplotAssay')
})

#' @rdname boxplotAssay
#' @export
setMethod('boxplotAssay', 'ProteinExperiment',
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
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }

  ## Data and options processing -----------------------------------------------
  mat <- assays(x)[[assayName]]

  loopCols <- .loopWrapper(x, 'conditionCol')

  ## make a data.frame for each condition
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      dfList <- list()
    }

    timeVec <- colData(x)[, giveMetaoption(x, 'timeCol')][loopCols[[i]]]

    values <- as.vector(mat[, loopCols[[i]]])
    tempDf <- data.frame(value = values,
                         time = rep(timeVec, each = nrow(mat)))
    tempDf$condition <- names(loopCols)[i]

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

#' @rdname boxplotAssay
#' @export
setMethod('boxplotAssay', 'PeptideExperiment',
          function(x,
                   assayName,
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  callNextMethod()

})

#' @rdname boxplotAssay
#' @export
setMethod('boxplotAssay', 'ProteomicsExperiment',
          function(x,
                   assayName,
                   mode = 'protein',
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  if (mode == 'protein') {

    boxplotAssay(x = x@ProteinExperiment,
              assayName = assayName,
              plotType = plotType,
              returnDataFrame = returnDataFrame,
              conditionCol = conditionCol,
              timeCol = timeCol)

  } else if (mode == 'peptide') {

    boxplotAssay(x = x@PeptideExperiment,
              assayName = assayName,
              plotType = plotType,
              returnDataFrame = returnDataFrame,
              conditionCol = conditionCol,
              timeCol = timeCol)
  }

})
