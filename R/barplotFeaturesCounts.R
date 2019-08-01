#' @rdname barplotFeaturesCounts
#' @name barplotFeaturesCounts
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
#' @return A ggplot2 barplot object or a \code{data.frame}.
#' @import ggplot2
#' @export
#' @examples
#' barplotFeaturesCounts(wormsPE, assayName = 'ratio')
setGeneric('barplotFeaturesCounts', function(x, ...){
  standardGeneric('barplotFeaturesCounts')
})


#' @rdname barplotFeaturesCounts
#' @export
setMethod('barplotFeaturesCounts',
          'ProteinExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
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
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  plotDf <- tryCatch(
    {
      colname <- giveMetaoption(x, 'conditionCol')
      colnames(plotDf)[colnames(plotDf) == colname] <- 'conditionCol'
      plotDf
    },
    error = function(cond){
      plotDf$conditionCol <- NA
      plotDf
    }
  )

  plotDf$rownames <- factor(rownames(plotDf), levels = rownames(plotDf))
  ## early return without plot
  if (returnDataFrame) {
    return(plotDf)
  }

  ## plot with fill depending if we have the conditionCol column
  if (all(is.na(plotDf$conditionCol))) {
    ggplot(data = plotDf,
           aes_string(x = 'rownames',
                      y = 'counts')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      theme_bw()
  } else {

    colname <- giveMetaoption(x, 'conditionCol')
    oldname <- colnames(plotDf)[colnames(plotDf) == colname]

    ggplot(data = plotDf,
           aes_string(x = 'rownames',
                      y = 'counts', fill = 'conditionCol')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)+
      labs(fill = oldname) +
      theme_bw()
  }

})



#' @rdname barplotFeaturesCounts
#' @export
setMethod('barplotFeaturesCounts',
          'PeptideExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  callNextMethod()

})

#' @rdname barplotFeaturesCounts
#' @export
setMethod('barplotFeaturesCounts',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  protPart <- barplotFeaturesCounts(x = x@ProteinExperiment,
                                    assayName = assayName,
                                    return = TRUE,
                                    conditionCol = conditionCol)

  peptPart <- barplotFeaturesCounts(x = x@PeptideExperiment,
                                    assayName = assayName,
                                    return = TRUE,
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
  if (returnDataFrame) {
    return(plotDf)
  }

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  if (all(is.na(plotDf$conditionCol))) {
    ggplot(data = plotDf,
           aes_string(x = 'sample',
                      y = 'counts')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  } else {
    ggplot(data = plotDf,
           aes_string(x = 'sample',
                      y = 'counts', fill = 'conditionCol')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      labs(fill = oldname) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  }


})
