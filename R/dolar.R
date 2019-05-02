
#' @export
setMethod('$', 'ProteomicsExperiment', function(x, name) {
  return(colData(x)[[name]])
})

#' @export
setReplaceMethod('$', 'ProteomicsExperiment', function(x, name, value) {
  colData(x)[[name]] <- value
  return(x)
})
