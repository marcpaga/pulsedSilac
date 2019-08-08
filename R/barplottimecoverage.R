#' @rdname barplotTimeCoverage
#' @name barplotTimeCoverage
#'
#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample. Anything
#' else than NA is considered detected.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param assayName Name of the assay to use in the plot.
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param ... Unused.
#'
#' @return A barplot or a \code{data.frame}.
#'
#' @examples
#'
#' barplotTimeCoverage(wormsPE, assayName = 'ratio')
#'
#' @import ggplot2
#' @export
setGeneric('barplotTimeCoverage', function(x, ...){
  standardGeneric('barplotTimeCoverage')
})

#' @rdname barplotTimeCoverage
#' @export
setMethod('barplotTimeCoverage',
          'ProteinExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  ## argument checking ---------------------------------------------------------

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## if metaoption given as argument put them in the object
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## data processing -----------------------------------------------------------

  ## assay used for the plotting
  mat <- assays(x)[[assayName]]

  ## get which columns belong to each condition
  loopCols <- .loopWrapper(x, 'conditionCol')

  ## count missing non-missing values per condition
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      plotDfList <- list()
    }
    plotDf <- data.frame()

    counts <- apply(mat[, loopCols[[i]]], 1, function(x) sum(!is.na(x)))
    plotDf <- data.frame(table(counts))
    if (is.null(names(loopCols))) {
      plotDf$condition <- as.factor(NA)
    } else {
      plotDf$condition <- as.factor(names(loopCols)[i])
    }
    plotDfList[[i]] <- plotDf

    if (i == length(loopCols)) {
      plotDf <- do.call('rbind', plotDfList)
    }
  }

  ## early return without plot
  if (returnDataFrame) {
    return(plotDf)
  }

  ## plotting ------------------------------------------------------------------

  ## there is one condition or not specified
  if (all(is.na(plotDf$condition))) {
    ggplot(data = plotDf,
           mapping = aes_string(x = 'counts',
                                y = 'Freq')) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      theme_bw()

  ## there is more than one condition
  } else {
    ggplot(data = plotDf,
           mapping = aes_string(x = 'counts',
                                y = 'Freq',
                                fill = 'condition')) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      theme_bw()
  }

})

#' @rdname barplotTimeCoverage
#' @export
setMethod('barplotTimeCoverage',
          'PeptideExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {


  callNextMethod()

})

#' @rdname barplotTimeCoverage
#' @export
setMethod('barplotTimeCoverage',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  protPart <- barplotTimeCoverage(x = x@ProteinExperiment,
                                  assayName = assayName,
                                  returnDataFrame = TRUE,
                                  conditionCol = conditionCol)

  peptPart <- barplotTimeCoverage(x = x@PeptideExperiment,
                                  assayName = assayName,
                                  returnDataFrame = TRUE,
                                  conditionCol = conditionCol)


  ## join the data.frames
  protPart$mode <- 'Protein'
  peptPart$mode <- 'Peptide'
  plotDf <- rbind(protPart, peptPart)
  plotDf$mode <- factor(plotDf$mode, levels = c('Protein', 'Peptide'))

  ## fix unordered counts when not all possibilities are there
  plotDf$counts <- factor(plotDf$counts,
                          levels = sort(unique(as.numeric(as.character(plotDf$counts)))))

  ## early return with no plot
  if (returnDataFrame) {
    return(plotDf)
  }

  ## only one condition or not specified
  if (all(is.na(plotDf$condition))) {

    ggplot(data = plotDf,
           mapping = aes_string(x = 'counts',
                                y = 'Freq')) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()

  ## more than one condition
  } else {
    ggplot(data = plotDf,
           mapping = aes_string(x = 'counts',
                                y = 'Freq',
                                fill = 'condition')) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  }


})
