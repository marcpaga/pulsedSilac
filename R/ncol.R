#' @rdname ProteomicsExperiment-accessors
#' @aliases ncol
#' @name ncol,ProteomicsExperiment-method
#' @export
#' @importFrom SummarizedExperiment colData
setMethod('ncol', 'ProteomicsExperiment', function(x){

  return(nrow(colData(x)))

})
