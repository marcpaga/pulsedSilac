#' @rdname SilacProteomicsExperiment-accessors
#' @aliases metadata,SilacProteomicsExperiment-method
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata', 'SilacProteomicsExperiment', function(x, ...){

  return(x@metadata)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases metadata<-,SilacProteomicsExperiment-method
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata<-', 'SilacProteomicsExperiment', function(x, ..., value){

  x@metadata <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})
