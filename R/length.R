#' @export
setMethod('length', 'ProteinExperiment', function(x){

  return(callNextMethod())

})


#' @export
setMethod('length', 'PeptideExperiment', function(x){

  return(callNextMethod())

})


#' @export
setMethod('length', 'ProteomicsExperiment', function(x){

  return(c(protein = length(x@ProteinExperiment),
           peptide = length(x@PeptideExperiment)))

})
