#' @rdname ProteomicsExperiment-accessors
#' @aliases nrow,ProteomicsExperiment-method
#' @usage NULL
#' @export
setMethod('nrow', 'ProteomicsExperiment', function(x){

  return(c(protein = nrow(x@ProteinExperiment),
           peptide = nrow(x@PeptideExperiment)))

})
