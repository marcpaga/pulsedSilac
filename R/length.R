#' @rdname SilacProteomicsExperiment-accessors
#' @aliases length,SilacProteomicsExperiment-method
#' @export
setMethod('length', 'SilacProteomicsExperiment', function(x){

  return(c(protein = length(x@SilacProteinExperiment),
           peptide = length(x@SilacPeptideExperiment)))

})
