#' @rdname ProteomicsExperiment-accessors
#' @aliases dim,ProteomicsExperiment-method
#' @export
setMethod('dim', 'ProteomicsExperiment', function(x){

  return(matrix(data = c(dim(x@ProteinExperiment), dim(x@PeptideExperiment)),
                ncol = 2,
                byrow = T,
                dimnames = list(c('protein', 'peptide'), c('row', 'column'))))

})
