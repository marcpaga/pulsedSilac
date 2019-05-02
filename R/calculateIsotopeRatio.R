#' @export
setGeneric('calculateIsotopeRatio', function(x, assayNameA, assayNameB, ...){
  standardGeneric('calculateIsotopeRatio')
})

#' @export
setMethod('calculateIsotopeRatio', 'ProteinExperiment',
          function(x, assayNameA = 'heavy_intensity',
                   assayNameB = 'light_intensity') {

  ratio_assay <- assays(x)[[assayNameA]]/assays(x)[[assayNameB]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeRatio', 'PeptideExperiment',
          function(x, assayNameA = 'heavy_intenisty',
                   assayNameB = 'light_intensity') {

  ratio_assay <- assays(x)[[assayNameA]]/assays(x)[[assayNameB]]
  assays(x)[['ratio']] <- ratio_assay
  return(x)

})

#' @export
setMethod('calculateIsotopeRatio', 'ProteomicsExperiment',
          function(x, assayNameA = 'heavy_intenisty',
                   assayNameB = 'light_intensity') {

  if (length(assayNameA) == 1) {
    assayNameA <- rep(assayNameA, 2)
  }

  if (length(assayNameB) == 1) {
    assayNameB <- rep(assayNameB, 2)
  }

  x@ProteinExperiment <- calculateIsotopeRatio(x@ProteinExperiment,
                                               assayNameA = assayNameA[1],
                                               assayNameB = assayNameB[1])

  x@PeptideExperiment <- calculateIsotopeRatio(x@PeptideExperiment,
                                               assayNameA = assayNameA[2],
                                               assayNameB = assayNameB[2])

  return(x)

})
