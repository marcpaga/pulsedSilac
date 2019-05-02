#' @export
setMethod('colnames', 'ProteomicsExperiment', function(x) {

  return(rownames(colData(x)))
})

#' @export
setMethod('colnames<-', 'ProteomicsExperiment', function(x, value) {

  rownames(colData(x)) <- value
  return(x)

})
