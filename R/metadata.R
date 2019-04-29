#' @title Metadata getter
#' @rdname metadata
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata', 'ProteomicsExperiment', function(x, ...){

  return(x@metadata)

})

#' @title Metadata setter
#' @rdname metadata
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata<-', 'ProteomicsExperiment', function(x, ..., value){

  x@metadata <- value
  validObject(x)
  return(x)

})
