#' @rdname classAccessors
#' @aliases dim-ProteinExperiment
#' @usage NULL
#' @export
setMethod('dim', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases dim-PeptideExperiment
#' @usage NULL
#' @export
setMethod('dim', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases dim-ProteomicsExperiment
#' @usage NULL
#' @export
setMethod('dim', 'ProteomicsExperiment', function(x){

  return(matrix(data = c(dim(x@ProteinExperiment), dim(x@PeptideExperiment)),
                ncol = 2,
                byrow = T,
                dimnames = list(c('protein', 'peptide'), c('row', 'column'))))

})
