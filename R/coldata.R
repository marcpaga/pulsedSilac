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

#' @export
setMethod('colData<-', 'PeptideExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @export
setMethod('colData', 'ProteomicsExperiment', function(x, ...) {

  return(x@colData)

})

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
