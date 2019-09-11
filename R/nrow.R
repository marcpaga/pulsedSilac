#' @rdname SilacProteomicsExperiment-accessors
#' @aliases nrow,SilacProteomicsExperiment-method
#' @usage NULL
#' @export
setMethod('nrow', 'SilacProteomicsExperiment', function(x){

  return(c(protein = nrow(x@SilacProteinExperiment),
           peptide = nrow(x@SilacPeptideExperiment)))

})
