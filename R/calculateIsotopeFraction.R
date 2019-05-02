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
setGeneric('calculateIsotopeFraction', function(x, assayName, ...){
  standardGeneric('calculateIsotopeFraction')
})

#' @export
setMethod('calculateIsotopeFraction', 'ProteinExperiment',
          function(x, assayName = 'ratio') {

  fraction_assay <- assays(x)[[assayName]]/
    ( 1 + assays(x)[[assayName]] )
  assays(x)[['fraction']] < fraction_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeFraction', 'PeptideExperiment',
          function(x, assayName = 'ratio') {

  fraction_assay <- assays(x)[[assayName]]/
    ( 1 + assays(x)[[assayName]] )
  assays(x)[['fraction']] < fraction_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeFraction', 'ProteomicsExperiment',
          function(x, assayName = 'ratio') {

  if (length(assayName) == 1) {
    assayName <- rep(assayName, 2)
  }

  x@ProteinExperiment <- calculateIsotopeFraction(x@ProteinExperiment,
                                                  assayName = assayName[1])

  x@PeptideExperiment <- calculateIsotopeFraction(x@PeptideExperiment,
                                                  assayName = assayName[2])

  return(x)

})
