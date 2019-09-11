#' @rdname SilacProteomicsExperiment-accessors
#' @aliases colnames,SilacProteomicsExperiment-method
#' @export
setMethod('colnames', 'SilacProteomicsExperiment', function(x) {

  return(rownames(x@colData))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases colnames<-,SilacProteomicsExperiment-method
#' @export
setMethod('colnames<-', 'SilacProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Colnames must be unique')
  }

  rownames(colData(x@SilacProteinExperiment)) <- value
  rownames(colData(x@SilacPeptideExperiment)) <- value
  validObject(x)
  return(x)

})
