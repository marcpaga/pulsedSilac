#' @rdname calculateOldIsotopePool
#' @name calculateOldIsotopePool
#' @title Estimates the fraction of old isotope for each time point
#'
#' To estimate how much of the "old" isotope is being used in "new" proteins we
#' can use the expression level of miss-cleaved peptides that contain a mix of
#' isotopes (one old and one new) and miss-cleaved peptides that contain only
#' new isotopes. This can be done using the following formula:
#' \deqn{\frac{1}{2\frac{Intensity_{lys8lys8}}{Intensity_{lys8lys0}} + 1} = lys_0 (Fraction)}
#' Which gives an idea of how much recyling (turnover understimation) is
#' happening.
#'
#' Both peptide types, mix of old/new isotope and two new isotopes, have to be
#' found in a time point to calculate the fraction of old isotope.
#'
#' @param x A \code{PeptideExperiment} or \code{ProteomicsExperiment} object.
#' @param newIsotopeAssayName \code{character} indicating the assay that
#' contains quantification data for miss-cleaved peptides with two new isotopes
#' incorporated.
#' @param mixIsotopeAssayName \code{character} indicating the assay that
#' contains quantification data for miss-cleaved peptides with one old isotope
#' and one new isotope incorporated.
#' @param ... Unused.
#'
#' @return A \code{PeptideExperiment} or \code{ProteomicsExperiment} with a
#' peptide assay entry named "oldIsotopePool".
#'
#' @examples
#' protPE <- ProtExp(wormsPE)
#' missPE <- addMisscleavedPeptides(x = protPE,
#'                                  newdata = recycleLightLysine,
#'                                  idColPept = 'Sequence',
#'                                  modCol = 'Modifications',
#'                                  dataCols = c(18:31))
#'
#' names(assays(missPE))[1:2] <- c('int_lys8lys8', 'int_lys8lys0')
#' missPE <- calculateOldIsotopePool(x = missPE, 'int_lys8lys8', 'int_lys8lys0')
#'
#' plotDistributionAssay(missPE, assayName = 'oldIsotopePool')
#'
#' @export
setGeneric('calculateOldIsotopePool', function(x, ...){
  standardGeneric('calculateOldIsotopePool')
})

#' @rdname calculateOldIsotopePool
#' @export
setMethod('calculateOldIsotopePool',
          'PeptideExperiment',
          function(x,
                   newIsotopeAssayName,
                   mixIsotopeAssayName) {

  ## argument check
  if (!(newIsotopeAssayName %in% assayNames(x)) |
      !(mixIsotopeAssayName %in% assayNames(x))) {
    stop('The given assayNames cannot be found in the assays')
  }

  newMat <- assays(x)[[newIsotopeAssayName]]
  mixMat <- assays(x)[[mixIsotopeAssayName]]

  oldMat <- 1/((2*newMat/mixMat) + 1)

  assays(x)[['oldIsotopePool']] <- oldMat
  return(x)

})

#' @rdname calculateOldIsotopePool
#' @export
setMethod('calculateOldIsotopePool',
          'ProteomicsExperiment',
          function(x,
                   newIsotopeAssayName,
                   mixIsotopeAssayName) {

  new.pe <- calculateOldIsotopePool(x = x@PeptideExperiment,
                                    newIsotopeAssayName = newIsotopeAssayName,
                                    mixIsotopeAssayName = mixIsotopeAssayName)

  x@PeptideExperiment <- new.pe

  return(x)

})

