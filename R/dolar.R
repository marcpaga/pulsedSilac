#' @rdname ProteomicsExperiment-accessors
#' @aliases $,ProteomicsExperiment-method
#' @export
setMethod('$', 'ProteomicsExperiment', function(x, name) {
  return(colData(x)[[name]])
})

#' @rdname ProteomicsExperiment-accessors
#' @aliases $<-,ProteomicsExperiment-method
#' @export
setReplaceMethod('$', 'ProteomicsExperiment', function(x, name, value) {
  colData(x)[[name]] <- value
  return(x)
})
