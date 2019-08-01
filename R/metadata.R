#' @rdname classAccessors
#' @aliases metadata
#' @usage NULL
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata', 'ProteomicsExperiment', function(x, ...){

  return(x@metadata)

})

#' @rdname classAccessors
#' @aliases metadata<-
#' @usage NULL
#' @importMethodsFrom S4Vectors metadata
#' @export
setMethod('metadata<-', 'ProteomicsExperiment', function(x, ..., value){

  x@metadata <- value
  validObject(x)
  return(x)

})
