#' @rdname SilacProteomicsExperiment-accessors
#' @name linkerDf
#' @aliases linkerDf,SilacProteomicsExperiment-method
#' @export
setMethod('linkerDf', 'SilacProteomicsExperiment', function(x){

  return(x@linkerDf)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name linkerDf<-
#' @aliases linkerDf<-,SilacProteomicsExperiment-method
#' @export
setMethod('linkerDf<-', 'SilacProteomicsExperiment', function(x, value){

  x@linkerDf <- value
  validObject(x)
  return(x)

})
