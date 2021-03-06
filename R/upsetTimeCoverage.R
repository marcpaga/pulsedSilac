#' @rdname upsetTimeCoverage
#' @name upsetTimeCoverage
#' @title Number of detected features per sample
#'
#' @description How many proteins/peptides are detected in each sample. Anything
#' else than NA is considered detected.
#'
#' @param x A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or a
#' \code{SilacProteomicsExperiment} object.
#' @param assayName Name of the assay to use in the plot.
#' @param maxMissing A \code{numerical} indicating how many timepoints can a
#' protein/peptide miss.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param returnList A \code{logical} indicating if the \code{list}
#' used for the plot should be returned instead.
#' @param ... Further arguments passed to \code{upset()}.
#'
#' @return A barplot or a \code{data.frame}.
#'
#' @examples
#' data('wormsPE')
#' upsetTimeCoverage(x = ProtExp(wormsPE),
#'                 assayName = 'ratio',
#'                 maxMissing = 2)
#'
#' @importFrom UpSetR fromList upset
#' @importFrom cowplot plot_grid
#' @export
setGeneric('upsetTimeCoverage', function(x, ...){
  standardGeneric('upsetTimeCoverage')
})

#' @rdname upsetTimeCoverage
#' @export
setMethod('upsetTimeCoverage',
          'SilacProteinExperiment',
          function(x,
                   assayName,
                   conditionCol,
                   maxMissing = 0,
                   returnList = FALSE,
                   ...) {

  ## argument checker ----------------------------------------------------------
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }
  if (!missing(conditionCol)) {
    metadata(x)[['conditionCol']] <- conditionCol
  }

  ## data and options processing -----------------------------------------------
  mat <- assays(x)[[assayName]]

  loopCols <- .loopWrapper(x, 'conditionCol')
  if (length(loopCols) < 2) {
    stop('There is only one condition, overlaps cannot be made')
  }

  ## count missing values in each of the loopable conditions
  for (i in seq_along(loopCols)) {
    if (i == 1) {
      logMatrix <- matrix(data = FALSE,
                          ncol = length(loopCols),
                          nrow = nrow(x))
    }

    nacounts <- apply(mat[, loopCols[[i]]], 1, function(x) sum(is.na(x)))
    ## TRUE if there are less or equal missing as maxMissing
    logMatrix[, i] <- (nacounts <= maxMissing)

  }

  ## get which are ones that pass the maxMissing
  overlapList <- apply(logMatrix, 2, which)
  if (is.matrix(overlapList)) {
    overlapList <- as.list(as.data.frame(overlapList))
  }

  names(overlapList) <- names(loopCols)

  ## no plot
  if (returnList) {
    return(overlapList)
  }

  upset(fromList(overlapList), ...)

})

#' @rdname upsetTimeCoverage
#' @export
setMethod('upsetTimeCoverage',
          'SilacPeptideExperiment',
          function(x,
                   assayName,
                   maxMissing = 0,
                   conditionCol,
                   returnList = FALSE,
                   ...) {

  callNextMethod()

})

#' @rdname upsetTimeCoverage
#' @importFrom grid grid.edit
#' @importFrom grid grid.grab
#' @export
setMethod('upsetTimeCoverage',
          'SilacProteomicsExperiment',
          function(x,
                   assayName,
                   maxMissing = 0,
                   conditionCol,
                   returnList = FALSE,
                   ...) {

  if (returnList) {
    prot_plot <- upsetTimeCoverage(x@SilacProteinExperiment,
                                 assayName = assayName,
                                 maxMissing = maxMissing,
                                 conditionCol = conditionCol,
                                 returnList = returnList,
                                 ... = ...)

    pept_plot <- upsetTimeCoverage(x@SilacPeptideExperiment,
                                 assayName = assayName,
                                 maxMissing = maxMissing,
                                 conditionCol = conditionCol,
                                 returnList = returnList,
                                 ... = ...)

    outList <- list(protein = prot_plot,
                    peptide = pept_plot)
    return(outList)
  }

  ## upsetPlots cannot be saved into a variable, therefore this grid trick
  ## has to be used for side by side plotting

  upsetTimeCoverage(x@SilacProteinExperiment,
                  assayName = assayName,
                  maxMissing = maxMissing,
                  conditionCol = conditionCol,
                  returnList = returnList,
                  ... = ...)
  grid.edit('arrange', name = 'arrange2', redraw = FALSE)
  vp1 <- grid.grab()

  upsetTimeCoverage(x@SilacPeptideExperiment,
                  assayName = assayName,
                  maxMissing = maxMissing,
                  conditionCol = conditionCol,
                  returnList = returnList,
                  ... = ...)

  grid.edit('arrange', name = 'arrange2', redraw = FALSE)
  vp2 <- grid.grab()


  plot_grid(vp1, vp2, labels = c('Protein', 'Peptide'), align = 'h')

})
