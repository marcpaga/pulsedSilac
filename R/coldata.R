#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases colData<-,SilacProteinExperiment-method
#' @export
#' @importFrom S4Vectors DataFrame
setMethod('colData<-', 'SilacProteinExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases colData<-,SilacPeptideExperiment-method
#' @export
#' @importFrom S4Vectors DataFrame
setMethod('colData<-', 'SilacPeptideExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }
  x@colData <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases colData,SilacProteomicsExperiment-method
#' @export
setMethod('colData', 'SilacProteomicsExperiment', function(x, ...) {

  return(x@colData)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases colData<-,SilacProteomicsExperiment-method
#' @export
setMethod('colData<-', 'SilacProteomicsExperiment', function(x, ..., value) {

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  x@colData <- value
  colData(x@SilacProteinExperiment) <- value
  colData(x@SilacPeptideExperiment) <- value
  validObject(x)
  return(x)

})
