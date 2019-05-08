#' @export
setGeneric('barplotTimeCoverage', function(x, ...){
  standardGeneric('barplotTimeCoverage')
})


#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample. Anything
#' else than NA is considered detected.
#'
#' @return A barplot or a \code{data.frame}.
#' @export
#' @import ggplot2
setMethod('barplotTimeCoverage',
          'ProteinExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol,
                   replicateTimeCol) {

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
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
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

      plotCol <- paste(colData(x)[, condCol],
                       colData(x)[, timeRepCol], sep = '.')

      plotCol <- insertElems(plotCol,
                             pos = vapply(loopCols, '[[', 1, 1),
                             elems = unique(plotCol))

      list(loopCols = loopCols, plotCol = plotCol)
    },
    error = function(c){
      tryCatch(
        {
          loopCols <- experimentLoopWrapper(x, 'cond')
          plotCol <- condCol

          plotCol <- insertElems(plotCol,
                                 pos = vapply(loopCols, '[[', 1, 1),
                                 elems = unique(plotCol))

          list(loopCols = loopCols, plotCol = plotCol)
        },
        error = function(c){
          loopCols <- list(seq_len(ncol(x)))
          plotCol <- rep(NA, ncol(x) + 1)
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

    plotDfList[[i]] <- plotDf

    if (i == length(loopCols)) {
      plotDf <- do.call('rbind', plotDfList)
      plotDf$group <- outList[[2]]
    }

  }


  ## early return without plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  ## cb palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$group))) {
    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq)) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)
  } else {

    colname <- giveMetaoption(x, 'conditionCol')
    oldname <- colnames(plotDf)[colnames(plotDf) == colname]

    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq, fill = group)) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname)
  }

})



#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample. Anything
#' else than NA is considered detected.
#'
#' @return A barplot or a \code{data.frame}.
#' @export
#' @import ggplot2
setMethod('barplotTimeCoverage',
          'PeptideExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol,
                   replicateTimeCol) {

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
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
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

      plotCol <- paste(colData(x)[, condCol],
                       colData(x)[, timeRepCol], sep = '.')

      plotCol <- insertElems(plotCol,
                             pos = vapply(loopCols, '[[', 1, 1),
                             elems = unique(plotCol))

      list(loopCols = loopCols, plotCol = plotCol)
    },
    error = function(c){
      tryCatch(
        {
          loopCols <- experimentLoopWrapper(x, 'cond')
          plotCol <- condCol

          plotCol <- insertElems(plotCol,
                                 pos = vapply(loopCols, '[[', 1, 1),
                                 elems = unique(plotCol))

          list(loopCols = loopCols, plotCol = plotCol)
        },
        error = function(c){
          loopCols <- list(seq_len(ncol(x)))
          plotCol <- rep(NA, ncol(x) + 1)
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

    plotDfList[[i]] <- plotDf

    if (i == length(loopCols)) {
      plotDf <- do.call('rbind', plotDfList)
      plotDf$group <- outList[[2]]
    }

  }

  ## early return without plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  ## cb palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$group))) {
    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq)) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)
  } else {

    colname <- giveMetaoption(x, 'conditionCol')
    oldname <- colnames(plotDf)[colnames(plotDf) == colname]

    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq, fill = group)) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname)
  }
})

#' @export
setMethod('barplotTimeCoverage',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol,
                   replicateTimeCol) {

  protPart <- barplotTimeCoverage(x = x@ProteinExperiment,
                                    assayName = assayName,
                                    return = 'data.frame',
                                    conditionCol = conditionCol,
                                    replicateTimeCol = replicateTimeCol)

  peptPart <- barplotTimeCoverage(x = x@PeptideExperiment,
                                    assayName = assayName,
                                    return = 'data.frame',
                                    conditionCol = conditionCol,
                                    replicateTimeCol = replicateTimeCol)

  ## join the data.frames
  protPart$mode <- 'Protein'
  peptPart$mode <- 'Peptide'
  plotDf <- rbind(protPart, peptPart)
  plotDf$mode <- factor(plotDf$mode, levels = c('Protein', 'Peptide'))

  ## early return with no plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  if (all(is.na(plotDf$group))) {
    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq)) +
      geom_bar(stat = 'identity') +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free')
  } else {
    ggplot(data = plotDf,
           aes(x = counts,
               y = Freq, fill = group)) +
      geom_bar(stat = 'identity', position = position_dodge()) +
      xlab('Timepoints') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free')
  }


})
