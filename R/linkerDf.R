#' @export
linkerDf <- setMethod('linkerDf', 'ProteomicsExperiment', function(x){

  return(x@linkerDf)

})


#' @export
linkerDf <- setMethod('linkerDf<-', 'ProteomicsExperiment', function(x){

  x@linkerDf <- value
  validObject(x)
  return(x)

})
