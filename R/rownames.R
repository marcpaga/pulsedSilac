#' @rdname ProteomicsExperiment-accessors
#' @name rownamesProt
#' @aliases rownamesProt,ProteomicsExperiment-method
#' @export
setMethod('rownamesProt', 'ProteomicsExperiment', function(x) {

  return(rownames(rowDataProt(x)))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rownamesProt<-
#' @aliases rownamesProt<-,ProteomicsExperiment-method
#' @export
setMethod('rownamesProt<-', 'ProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Rownames must be unique')
  }

  rownames(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name rownamesPept
#' @aliases rownamesPept,ProteomicsExperiment-method
#' @export
setMethod('rownamesPept', 'ProteomicsExperiment', function(x) {

  return(rownames(rowDataPept(x)))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rownamesPept<-
#' @aliases rownamesPept<-,ProteomicsExperiment-method
#' @export
setMethod('rownamesPept<-', 'ProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Rownames must be unique')
  }

  rownames(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
