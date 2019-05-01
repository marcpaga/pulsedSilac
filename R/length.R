#' @export
setMethod('length', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})


#' @export
setMethod('length', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})


#' @export
setMethod('length', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(c(protein = length(x@ProteinExperiment),
           peptide = length(x@PeptideExperiment)))

})
