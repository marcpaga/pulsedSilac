#' @export
setGeneric('plotAssay', function(x, ...){
  standardGeneric('plotAssay')
})

#' @title Distribution of assay data per condition and timepoint.
#'
#' @description Plot the distribution of the data stored in
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
#' @param replicateTimeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different time replicates.
#'
#'
#' @importFrom ggridges geom_density_ridges
#' @import ggplot2
#' @export
setMethod('plotAssay', 'ProteinExperiment',
          function(x,
                   assayName,
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol,
                   replicateTimeCol) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  if (!plotType %in% c('boxplot', 'density')) {
    txt <- c('plotType must be "boxplot" or "density"')
    stop(txt)
  }

  ## count how many proteins per sample
  mat <- assays(x)[[assayName]]

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  ## first it tries to get both condition and time replicates, if it fails then
  ## only condition and if it fails then all the samples are grouped together
  outList <- tryCatch(
    {
      loopCols <- experimentLoopWrapper(x, 'cond.timerep')

      condCol <- giveMetaoption(x, 'conditionCol')
      timeRepCol <- giveMetaoption(x, 'replicateTimeCol')

      plotCol <- unique(paste(colData(x)[, condCol],
                              colData(x)[, timeRepCol], sep = '.'))

      list(loopCols = loopCols, plotCol = plotCol)
    },
    error = function(c){
      tryCatch(
        {
          loopCols <- experimentLoopWrapper(x, 'cond')
          condCol <- giveMetaoption(x, 'conditionCol')
          plotCol <- unique(colData(x)[, condCol])

          list(loopCols = loopCols, plotCol = plotCol)
        },
        error = function(c){
          'There is only 1 condition (?)'
        }
      )
    }
  )

  loopCols <- outList[[1]]


  for (i in seq_along(loopCols)) {
    if (i == 1) {
      dfList <- list()
    }

    timeVec <- colData(x)[,giveMetaoption(x, 'timeCol')][loopCols[[i]]]

    values <- as.vector(mat[,loopCols[[i]]])
    tempDf <- data.frame(value = values,
                         time = rep(timeVec, each = nrow(mat)))
    tempDf$condition <- outList$plotCol[i]

    dfList[[i]] <- tempDf
  }

  plotDf <- do.call('rbind', dfList)
  plotDf$time <- as.factor(plotDf$time)
  plotDf$condition <- as.factor(plotDf$condition)

  if (returnDataFrame) {
    colnames(plotDf)[1] <- assayName
    return(plotDf)
  }

  if (plotType == 'density') {

    p <- ggplot(data = plotDf) +
      geom_density_ridges(aes_string(x = 'value',
                                     y = 'time',
                                     fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~condition) +
      labs(x = assayName) +
      theme_classic()

  } else if (plotType == 'boxplot') {

    p <- ggplot(data = plotDf) +
      geom_boxplot(aes_string(x = 'time',
                              y = 'value',
                              fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      labs(y = assayName) +
      theme_classic()

  }

  p

})

#' @export
setMethod('plotAssay', 'PeptideExperiment',
          function(x,
                   assayName,
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol,
                   replicateTimeCol) {

  callNextMethod()

})

#' @export
setMethod('plotAssay', 'ProteomicsExperiment',
          function(x,
                   assayName,
                   mode = 'protein',
                   plotType = 'boxplot',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol,
                   replicateTimeCol) {

  if (mode == 'protein') {

    plotAssay(x = x@ProteinExperiment,
              assayName = assayName,
              plotType = plotType,
              returnDataFrame = returnDataFrame,
              conditionCol = conditionCol,
              timeCol = timeCol,
              replicateTimeCol = replicateTimeCol)

  } else if (mode == 'peptide') {

    plotAssay(x = x@PeptideExperiment,
              assayName = assayName,
              plotType = plotType,
              returnDataFrame = returnDataFrame,
              conditionCol = conditionCol,
              timeCol = timeCol,
              replicateTimeCol = replicateTimeCol)
  }

})
