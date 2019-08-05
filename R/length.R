#' @rdname ProteomicsExperiment-accessors
#' @aliases length,ProteomicsExperiment-method
#' @export
setMethod('length', 'ProteomicsExperiment', function(x){

  return(c(protein = length(x@ProteinExperiment),
           peptide = length(x@PeptideExperiment)))

})
