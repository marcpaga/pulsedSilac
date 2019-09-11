#' @rdname SilacProteomicsExperiment-accessors
#' @name rownamesProt
#' @aliases rownamesProt,SilacProteomicsExperiment-method
#' @export
setMethod('rownamesProt', 'SilacProteomicsExperiment', function(x) {

  return(rownames(rowDataProt(x)))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rownamesProt<-
#' @aliases rownamesProt<-,SilacProteomicsExperiment-method
#' @export
setMethod('rownamesProt<-', 'SilacProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Rownames must be unique')
  }

  rownames(x@SilacProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rownamesPept
#' @aliases rownamesPept,SilacProteomicsExperiment-method
#' @export
setMethod('rownamesPept', 'SilacProteomicsExperiment', function(x) {

  return(rownames(rowDataPept(x)))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rownamesPept<-
#' @aliases rownamesPept<-,SilacProteomicsExperiment-method
#' @export
setMethod('rownamesPept<-', 'SilacProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Rownames must be unique')
  }

  rownames(x@SilacPeptideExperiment) <- value
  validObject(x)
  return(x)

})
