#' @rdname barplotTimeCoverage
#' @name barplotTimeCoverage
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

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## assay used for the plotting
  mat <- assays(x)[[assayName]]

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
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
          loopCols <- list(seq_len(ncol(x)))
          plotCol <- NA
          list(loopCols = loopCols, plotCol = plotCol)
        }
      )
    }
  )

  loopCols <- outList[[1]]

  for (i in seq_along(loopCols)) {
    if (i == 1) {
      plotDfList <- list()
    }

    plotDf <- data.frame()

    counts <- apply(mat[, loopCols[[i]]], 1, function(x) sum(!is.na(x)))
    plotDf <- data.frame(table(counts))
    plotDf$group <- rep(outList[[2]][i], nrow(plotDf))

    plotDfList[[i]] <- plotDf

    if (i == length(loopCols)) {
      plotDf <- do.call('rbind', plotDfList)
    }

  }


  ## early return without plot
  if (returnDataFrame) {
    return(plotDf)
  }

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$group))) {
    ggplot(data = plotDf,
           aes_string(x = 'counts',
                      y = 'Freq')) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      theme_bw()
  } else {

    colname <- giveMetaoption(x, 'conditionCol')
    oldname <- colnames(plotDf)[colnames(plotDf) == colname]

    ggplot(data = plotDf,
           aes_string(x = 'counts',
                     y = 'Freq', fill = 'group')) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname) +
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

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (all(is.na(plotDf$group))) {
    ggplot(data = plotDf,
           aes_string(x = 'counts',
                      y = 'Freq')) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  } else {
    ggplot(data = plotDf,
           aes_string(x = 'counts',
                      y = 'Freq', fill = 'group')) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  }


})
