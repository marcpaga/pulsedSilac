#' @rdname calculateIsotopeRatio
#' @name calculateIsotopeRatio
#'
#' @title Ratio calculation
#'
#' @description Calculates the ratio between the new isotope and the old
#' isotope (new/old).
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param newIsotopeAssay A \code{character} indicating the assay name that
#' has the new isotope intensity data.
#' @param oldIsotopeAssay A \code{character} indicating the assay name that
#' has the old isotope intensity data.
#' @param ... Unused.
#'
#' @return A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object with an added assay named 'ratio'.
#'
#' @examples
#' data('wormsPE')
#' calculateIsotopeRatio(x = wormsPE,
#'                       newIsotopeAssay = 'int_heavy',
#'                       oldIsotopeAssay = 'int_light')
#'
#' @export
setGeneric('calculateIsotopeRatio', function(x, newIsotopeAssay,
                                             oldIsotopeAssay, ...){
  standardGeneric('calculateIsotopeRatio')
})

#' @rdname calculateIsotopeRatio
#' @export
setMethod('calculateIsotopeRatio', 'ProteinExperiment',
          function(x, newIsotopeAssay, oldIsotopeAssay) {

  ratio_assay <- assays(x)[[newIsotopeAssay]]/assays(x)[[oldIsotopeAssay]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @rdname calculateIsotopeRatio
#' @export
setMethod('calculateIsotopeRatio', 'PeptideExperiment',
          function(x, newIsotopeAssay, oldIsotopeAssay) {

  ratio_assay <- assays(x)[[newIsotopeAssay]]/assays(x)[[oldIsotopeAssay]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @rdname calculateIsotopeRatio
#' @export
setMethod('calculateIsotopeRatio', 'ProteomicsExperiment',
          function(x, newIsotopeAssay, oldIsotopeAssay) {

  if (length(newIsotopeAssay) == 1) {
    newIsotopeAssay <- rep(newIsotopeAssay, 2)
  }

  if (length(oldIsotopeAssay) == 1) {
    oldIsotopeAssay <- rep(oldIsotopeAssay, 2)
  }

  x@ProteinExperiment <- calculateIsotopeRatio(x@ProteinExperiment,
                                        newIsotopeAssay = newIsotopeAssay[1],
                                        oldIsotopeAssay = oldIsotopeAssay[1])

  x@PeptideExperiment <- calculateIsotopeRatio(x@PeptideExperiment,
                                        newIsotopeAssay = newIsotopeAssay[2],
                                        oldIsotopeAssay = oldIsotopeAssay[2])

  return(x)

})

