#' @rdname ProteomicsExperiment-accessors
#' @aliases metadata,ProteomicsExperiment-method
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata', 'ProteomicsExperiment', function(x, ...){

  return(x@metadata)

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases metadata<-,ProteomicsExperiment-method
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata<-', 'ProteomicsExperiment', function(x, ..., value){

  x@metadata <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})
