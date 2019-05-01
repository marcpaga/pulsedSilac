#' @export
setMethod('colData', 'ProteomicsExperiment', function(x, ...) {

  return(x@colData)

})

#' @export
setMethod('colData<-', 'ProteomicsExperiment', function(x, ..., value) {

  x@colData <- value
  validObject(x)
  return(x)

})
