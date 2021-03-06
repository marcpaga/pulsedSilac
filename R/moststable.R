#' @rdname mostStable
#' @name mostStable
#'
#' @title Most stable proteins/peptides
#'
#' @description Finds which are the most stable proteins/peptides across
#' the entire experiment. These proteins/peptides can be used to estimate
#' the cell growth of each condition.
#'
#' @details Proteins/peptides that are not found in all timepoints and in all
#' conditions are not considered. The stability is based on ranking
#' heavy label incorporation for each timepoint; therefore, lower values
#' are correlated to higher stability.
#'
#' @param x A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or a
#' \code{SilacProteomicsExperiment} object.
#' @param assayName Name of the assay to use.
#' @param n A \code{numeric} indicating how many proteins should be returned.
#' @param mode A \code{character} indicating which level of data to use,
#' either "protein" or "peptide". Only relevant for ProteomicsExperiment
#' inputs.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param ... Unused.
#'
#' @return A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or a
#' \code{SilacProteomicsExperiment} object with the n most stable
#' proteins/peptides.
#'
#' @examples
#' data('mefPE')
#' mostStable(mefPE, assayName = 'fraction', n = 50)
#'
#' @export
setGeneric('mostStable', function(x, ...){
  standardGeneric('mostStable')
})


#' @rdname mostStable
#' @export
setMethod('mostStable',
          'SilacProteinExperiment',
          function(x,
                   assayName,
                   n,
                   conditionCol) {


  ## argument checker ----------------------------------------------------------
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }
  if (!missing(conditionCol)) {
    metadata(x)[['conditionCol']] <- conditionCol
  }
  if (missing(n)) {
    stop('Must provide an "n" value')
  }

  ## data processing -----------------------------------------------------------
  mat <- assays(x)[[assayName]]

  loopCols <- .loopWrapper(x, 'conditionCol')

  ## initialize matrix were the stability ranks will go
  rankRes <- matrix(data = NA,
                    ncol = length(loopCols),
                    nrow = nrow(x))

  ## rank for each condition
  for (j in seq_along(loopCols)) {
    tempMat <- assays(x)[[assayName]][, loopCols[[j]]]
    rankMat <- apply(tempMat, 2, rank, na.last = 'keep')
    resMat <- rank(apply(rankMat, 1, sum, na.rm = FALSE), na.last = 'keep')
    rankRes[, j] <- resMat
  }

  ## apply a global rank across conditions and count the n top
  cutoff <- sort(apply(rankRes, 1, mean))[n]
  ## return the most stable ones
  stablePE <- x[which(apply(rankRes, 1, mean) <= cutoff), ]

  return(stablePE)

})

#' @rdname mostStable
#' @export
setMethod('mostStable',
          'SilacPeptideExperiment',
          function(x,
                   assayName,
                   n,
                   conditionCol) {

  callNextMethod()

})

#' @rdname mostStable
#' @export
setMethod('mostStable',
          'SilacProteomicsExperiment',
          function(x,
                   assayName,
                   n,
                   mode,
                   conditionCol) {

  if (mode == 'protein') {
    mostStable(x = x@SilacProteinExperiment,
               assayName = assayName,
               n = n,
               conditionCol = conditionCol)
  } else if (mode == 'peptide') {
    mostStable(x = x@SilacPeptideExperiment,
               assayName = assayName,
               n = n,
               conditionCol = conditionCol)
  }

})
