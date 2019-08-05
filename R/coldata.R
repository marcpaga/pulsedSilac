#' @rdname ProteinPeptideExperiment-accessors
#' @aliases colData<-,ProteinExperiment-method
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

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases colData<-,PeptideExperiment-method
#' @export
#' @importFrom S4Vectors DataFrame
setMethod('colData<-', 'PeptideExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases colData,ProteomicsExperiment-method
#' @export
setMethod('colData', 'ProteomicsExperiment', function(x, ...) {

  return(x@colData)

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases colData<-,ProteomicsExperiment-method
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
