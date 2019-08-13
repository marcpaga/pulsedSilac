#' @rdname barplotCounts
#' @name barplotCounts
#'
#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample.
#' \code{NA} are considered missing values.
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
#' @return A ggplot2 barplot object or a \code{data.frame}.
#'
#' @examples
#' data('wormsPE')
#' barplotCounts(wormsPE, assayName = 'ratio')
#'
#' @import ggplot2
#' @export
setGeneric('barplotCounts', function(x, ...){
  standardGeneric('barplotCounts')
})


#' @rdname barplotCounts
#' @export
setMethod('barplotCounts',
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

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## if metaoption given as argument put them in the object
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }

  ## data processing -----------------------------------------------------------

  ## count how many proteins per sample, NAs are missing values
  mat <- assays(x)[[assayName]]
  counts <- apply(mat, 2, function(x) sum(!is.na(x)))

  ## use colData as skeleton
  plotDf <- as.data.frame(colData(x))
  plotDf$counts <- counts

  ## process if there are multiple conditions and rename columns for plotting
  mOptColName <- .giveMetaoption(x, 'conditionCol')
  plotDf <- .processColDataPlot(plotDf, mOptColName, 'condition')
  plotDf$condition <- as.factor(plotDf$condition)

  plotDf$rownames <- factor(rownames(plotDf), levels = rownames(plotDf))
  ## early return without plot
  if (returnDataFrame) {
    return(plotDf)
  }

  ## plotting ------------------------------------------------------------------

  ## there is only one condition or is not provided
  if (all(is.na(plotDf[, 'condition']))) {
    ggplot(data = plotDf,
           mapping = aes_string(x = 'rownames',
                                y = 'counts')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)+
      theme_bw()
  ## there is more than one condition
  } else {

    ggplot(data = plotDf,
           mapping = aes_string(x = 'rownames',
                                y = 'counts',
                                fill = 'condition')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette)+
      theme_bw()
  }

})


#' @rdname barplotCounts
#' @export
setMethod('barplotCounts',
          'PeptideExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  callNextMethod()

})

#' @rdname barplotCounts
#' @export
setMethod('barplotCounts',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol) {

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  ## argument checking is done in the Protein and Peptide experiment functions

  protPart <- barplotCounts(x = x@ProteinExperiment,
                            assayName = assayName,
                            return = TRUE,
                            conditionCol = conditionCol)

  peptPart <- barplotCounts(x = x@PeptideExperiment,
                            assayName = assayName,
                            return = TRUE,
                            conditionCol = conditionCol)

  ## join the data.frames
  protPart$mode <- 'Protein'
  peptPart$mode <- 'Peptide'
  protPart$sample <- factor(rownames(protPart), levels = rownames(protPart))
  peptPart$sample <- factor(rownames(peptPart), levels = rownames(peptPart))
  plotDf <- rbind(protPart, peptPart)
  ## to do facet wrap
  plotDf$mode <- factor(plotDf$mode, levels = c('Protein', 'Peptide'))

  ## early return with no plot
  if (returnDataFrame) {
    return(plotDf)
  }

  if (all(is.na(plotDf$condition))) {
    ggplot(data = plotDf,
           mapping = aes_string(x = 'sample',
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
           mapping = aes_string(x = 'sample',
                                y = 'counts',
                                fill = 'condition')) +
      geom_bar(stat = 'identity') +
      xlab('Sample') +
      ylab('Counts') +
      theme(panel.border = element_rect(fill = NA)) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~mode, scales = 'free') +
      theme_bw()
  }


})
