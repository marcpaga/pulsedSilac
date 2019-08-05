#' @rdname ProteomicsExperiment-accessors
#' @name linkerDf
#' @aliases linkerDf,ProteomicsExperiment-method
#' @export
setMethod('linkerDf', 'ProteomicsExperiment', function(x){

  return(x@linkerDf)

})

#' @rdname ProteomicsExperiment-accessors
#' @name linkerDf<-
#' @aliases linkerDf<-,ProteomicsExperiment-method
#' @export
setMethod('linkerDf<-', 'ProteomicsExperiment', function(x, value){

  x@linkerDf <- value
  validObject(x)
  return(x)

})
