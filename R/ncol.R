#' @rdname classAccessors
#' @aliases ncol
#' @usage NULL
#' @export
setMethod('ncol', 'ProteinExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases ncol
#' @usage NULL
#' @export
setMethod('ncol', 'PeptideExperiment', function(x){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases ncol
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment colData
setMethod('ncol', 'ProteomicsExperiment', function(x){

  return(nrow(colData(x)))

})
