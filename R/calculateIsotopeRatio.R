#' @export
setGeneric('calculateIsotopeRatio', function(x, newIsotopeAssay, oldIsotopeAssay, ...){
  standardGeneric('calculateIsotopeRatio')
})

#' @export
setMethod('calculateIsotopeRatio', 'ProteinExperiment',
          function(x, newIsotopeAssay = 'heavy_intensity',
                   oldIsotopeAssay = 'light_intensity') {

  ratio_assay <- assays(x)[[newIsotopeAssay]]/assays(x)[[oldIsotopeAssay]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeRatio', 'PeptideExperiment',
          function(x, newIsotopeAssay = 'heavy_intenisty',
                   oldIsotopeAssay = 'light_intensity') {

  ratio_assay <- assays(x)[[newIsotopeAssay]]/assays(x)[[oldIsotopeAssay]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeRatio', 'ProteomicsExperiment',
          function(x, newIsotopeAssay = 'heavy_intenisty',
                   oldIsotopeAssay = 'light_intensity') {

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
