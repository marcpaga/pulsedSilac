#' @rdname SilacProteomicsExperiment-accessors
#' @aliases $,SilacProteomicsExperiment-method
#' @export
setMethod('$', 'SilacProteomicsExperiment', function(x, name) {
  return(colData(x)[[name]])
})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases $<-,SilacProteomicsExperiment-method
#' @export
setReplaceMethod('$', 'SilacProteomicsExperiment', function(x, name, value) {
  colData(x)[[name]] <- value
  return(x)
})
