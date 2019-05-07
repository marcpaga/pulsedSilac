#' @export
setGeneric('barplotFeaturesCounts', function(x, ...){
  standardGeneric('barplotFeaturesCounts')
})


#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample. Anything
#' else than NA is considered detected.
#'
#' @return A barplot or a \code{data.frame}.
#' @export
#' @import ggplot2
setMethod('barplotFeaturesCounts',
          'ProteinExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## count how many proteins per sample
  mat <- assays(x)[[assayName]]
  counts <- apply(mat, 2, function(x) sum(!is.na(x)))

  plotDf <- as.data.frame(colData(x))
  plotDf$counts <- counts

  ## cb palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  plotDfList <- tryCatch(
    {
      colname <- giveMetaoption(x, 'conditionCol')
      oldname <- colnames(plotDf)[colnames(plotDf) == colname]
      colnames(plotDf)[colnames(plotDf) == colname] <- 'conditionCol'
      list(plotDf, oldname)
    },
    error = function(cond){
      oldname <- NA
      plotDf$conditionCol <- NA
      list(plotDf, oldname)
    }
  )
  plotDf <- plotDfList[[1]]
  ## name in the legend
  oldname <- plotDfList[[2]]

  ## early return without plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$conditionCol))) {
    ggplot(data = plotDf,
           aes(x = factor(rownames(plotDf), levels = rownames(plotDf)),
               y = counts)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)
  } else {
    ggplot(data = plotDf,
           aes(x = factor(rownames(plotDf), levels = rownames(plotDf)),
               y = counts, fill = conditionCol)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)+
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
setMethod('barplotFeaturesCounts',
          'PeptideExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## count how many proteins per sample
  mat <- assays(x)[[assayName]]
  counts <- apply(mat, 2, function(x) sum(!is.na(x)))

  plotDf <- as.data.frame(colData(x))
  plotDf$counts <- counts

  ## cb palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  plotDfList <- tryCatch(
    {
      colname <- giveMetaoption(x, 'conditionCol')
      oldname <- colnames(plotDf)[colnames(plotDf) == colname]
      colnames(plotDf)[colnames(plotDf) == colname] <- 'conditionCol'
      list(plotDf, oldname)
    },
    error = function(cond){
      oldname <- NA
      plotDf$conditionCol <- NA
      list(plotDf, oldname)
    }
  )
  plotDf <- plotDfList[[1]]
  ## name in the legend
  oldname <- plotDfList[[2]]

  ## early return with no plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$conditionCol))) {
    ggplot(data = plotDf,
           aes(x = factor(rownames(plotDf), levels = rownames(plotDf)),
               y = counts)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)
  } else {
    ggplot(data = plotDf,
           aes(x = factor(rownames(plotDf), levels = rownames(plotDf)),
               y = counts, fill = conditionCol)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname)
  }

})

#' @export
setMethod('barplotFeaturesCounts',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   return = 'plot',
                   conditionCol) {

  protPart <- barplotFeaturesCounts(x = x@ProteinExperiment,
                                    assayName = assayName,
                                    return = 'data.frame',
                                    conditionCol = conditionCol)

  peptPart <- barplotFeaturesCounts(x = x@PeptideExperiment,
                                    assayName = assayName,
                                    return = 'data.frame',
                                    conditionCol = conditionCol)

  ## join the data.frames
  protPart$mode <- 'Protein'
  peptPart$mode <- 'Peptide'
  protPart$sample <- factor(rownames(protPart), levels = rownames(protPart))
  peptPart$sample <- factor(rownames(peptPart), levels = rownames(peptPart))
  plotDf <- rbind(protPart, peptPart)
  plotDf$mode <- factor(plotDf$mode, levels = c('Protein', 'Peptide'))

  ## name in the legend
  oldname <- colnames(colData(x))[which(colnames(plotDf) == 'conditionCol')]

  ## early return with no plot
  if (return == 'data.frame') {
    return(plotDf)
  }

  if (all(is.na(plotDf$conditionCol))) {
    ggplot(data = plotDf,
           aes(x = sample,
               y = counts)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free')
  } else {
    ggplot(data = plotDf,
           aes(x = sample,
               y = counts, fill = conditionCol)) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname) +
      facet_wrap(~mode, scales = 'free')
  }


})
