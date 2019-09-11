#' @rdname SilacProteomicsExperiment-accessors
#' @aliases dim,SilacProteomicsExperiment-method
#' @export
setMethod('dim', 'SilacProteomicsExperiment', function(x){

  return(matrix(data = c(dim(x@SilacProteinExperiment),
                         dim(x@SilacPeptideExperiment)),
                ncol = 2,
                byrow = T,
                dimnames = list(c('protein', 'peptide'), c('row', 'column'))))

})
