#' Calculates the incorporated isotope fraction
#'
#' Calculates the fraction of an isotope ratio using the followin formula:
#'
#' \deqn{Isotope fraction_{A} = \frac{ratio}{ratio + 1}}
#'
#' The ratio should be calculated as:
#' \deqn{ratio = isotope_{A}/isotope_{B}}
#'
#' @param x a \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object.
#' @param assayName a \code{character} with the name of the assay that contains
#' isotope ratios.
#'
#' @return a \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object with additional assays named "fraction".
#' @examples
#'
#' @seealso \code{\link{ProteomicsExperiment-class}}
#' @export
setGeneric('calculateIsotopeFraction', function(x, ratioAssay, ...){
  standardGeneric('calculateIsotopeFraction')
})

#' @export
setMethod('calculateIsotopeFraction', 'ProteinExperiment',
          function(x, ratioAssay = 'ratio') {

  fraction_assay <- assays(x)[[ratioAssay]]/( 1 + assays(x)[[ratioAssay]] )
  assays(x)[['fraction']] <- fraction_assay

  return(x)

})

#' @export
setMethod('calculateIsotopeFraction', 'PeptideExperiment',
          function(x, ratioAssay = 'ratio') {

  fraction_assay <- assays(x)[[ratioAssay]]/( 1 + assays(x)[[ratioAssay]] )
  assays(x)[['fraction']] <- fraction_assay

  return(x)

})

#' @export
setMethod('calculateIsotopeFraction', 'ProteomicsExperiment',
          function(x, ratioAssay = 'ratio') {

  if (length(ratioAssay) == 1) {
    ratioAssay <- rep(ratioAssay, 2)
  }

  x@ProteinExperiment <- calculateIsotopeFraction(x@ProteinExperiment,
                                                  ratioAssay = ratioAssay[1])

  x@PeptideExperiment <- calculateIsotopeFraction(x@PeptideExperiment,
                                                  ratioAssay = ratioAssay[2])

  return(x)

})
