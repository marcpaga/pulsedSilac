#' @rdname classAccessors
#' @aliases nrow
#' @usage NULL
#' @export
setMethod('nrow', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases nrow
#' @usage NULL
#' @export
setMethod('nrow', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases nrow
#' @usage NULL
#' @export
setMethod('nrow', 'ProteomicsExperiment', function(x){

  return(c(protein = nrow(x@ProteinExperiment),
           peptide = nrow(x@PeptideExperiment)))

})
