#' @rdname classAccessors
#' @usage NULL
#' @export
setMethod('linkerDf', 'ProteomicsExperiment', function(x){

  return(x@linkerDf)

})

#' @rdname classAccessors
#' @usage NULL
#' @export
setMethod('linkerDf<-', 'ProteomicsExperiment', function(x, value){

  x@linkerDf <- value
  validObject(x)
  return(x)

})
