#' @rdname SilacProteomicsExperiment-accessors
#' @aliases ncol
#' @name ncol,SilacProteomicsExperiment-method
#' @export
#' @importFrom SummarizedExperiment colData
setMethod('ncol', 'SilacProteomicsExperiment', function(x){

  return(nrow(colData(x)))

})
