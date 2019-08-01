#' @rdname classAccessors
#' @aliases colData<-
#' @usage NULL
#' @export
#' @importFrom S4Vectors DataFrame
setMethod('colData<-', 'ProteinExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @aliases colData<-
#' @usage NULL
#' @export
setMethod('colData<-', 'PeptideExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @aliases colData
#' @usage NULL
#' @export
setMethod('colData', 'ProteomicsExperiment', function(x, ...) {

  return(x@colData)

})

#' @rdname classAccessors
#' @aliases colData<-
#' @usage NULL
#' @export
setMethod('colData<-', 'ProteomicsExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  x@colData <- value
  colData(x@ProteinExperiment) <- value
  colData(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
