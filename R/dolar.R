#' @rdname classAccessors
#' @aliases $
#' @usage NULL
#' @export
setMethod('$', 'ProteomicsExperiment', function(x, name) {
  return(colData(x)[[name]])
})

#' @rdname classAccessors
#' @aliases $<-
#' @usage NULL
#' @export
setReplaceMethod('$', 'ProteomicsExperiment', function(x, name, value) {
  colData(x)[[name]] <- value
  return(x)
})
