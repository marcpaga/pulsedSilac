#' @export
setMethod('ncol', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('ncol', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @export
setMethod('ncol', 'ProteomicsExperiment', function(x){

  return(nrow(colData(x)))

})
