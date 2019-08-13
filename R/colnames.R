#' @rdname ProteomicsExperiment-accessors
#' @aliases colnames,ProteomicsExperiment-method
#' @export
setMethod('colnames', 'ProteomicsExperiment', function(x) {

  return(rownames(x@colData))

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases colnames<-,ProteomicsExperiment-method
#' @export
setMethod('colnames<-', 'ProteomicsExperiment', function(x, value) {

  if (length(value) != length(unique(value))) {
    stop('Colnames must be unique')
  }

  rownames(colData(x@ProteinExperiment)) <- value
  rownames(colData(x@PeptideExperiment)) <- value
  validObject(x)
  return(x)

})
