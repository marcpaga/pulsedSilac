#' @export
setMethod('linkerDf', 'ProteomicsExperiment', function(x){

  return(x@linkerDf)

})


#' @export
setMethod('linkerDf<-', 'ProteomicsExperiment', function(x){

  x@linkerDf <- value
  validObject(x)
  return(x)

})
