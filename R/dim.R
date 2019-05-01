#' @export
setMethod('dim', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('dim', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('dim', 'ProteomicsExperiment', function(x){

  return(matrix(c(dim(x@ProteinExperiment), dim(x@PeptideExperiment)),
                ncol = 2, byrow = T,
                dimnames = list(c('protein', 'peptide'), c('row', 'column'))))

})
