#' @export
setMethod('nrow', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('nrow', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('nrow', 'ProteomicsExperiment', function(x){

  return(c(protein = nrow(x@ProteinExperiment),
    peptide = nrow(x@PeptideExperiment)))

})
