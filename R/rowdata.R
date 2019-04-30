#' @export
setMethod('rowData', 'ProteinExperiment', function(x, use.names, ...){

  return(callNextMethod())

})

#' @export
setMethod('rowData<-', 'ProteinExperiment', function(x, ..., value){

  return(callNextMethod())

})

#' @export
setMethod('rowData', 'PeptideExperiment', function(x, use.names, ...){

  return(callNextMethod())

})

#' @export
setMethod('rowData<-', 'PeptideExperiment', function(x, ..., value){

  return(callNextMethod())

})

#' @export
setMethod('rowData', 'ProteomicsExperiment', function(x, use.names, ...){

  return(list(protein = rowData(x@ProteinExperiment),
              peptide = rowData(x@PeptideExperiment)))

})

#' @export
setMethod('rowData<-', 'ProteomicsExperiment', function(x, ..., value){

  if (!is.list(value) | length(value) != 2) {
    stop('Value must be a list of length 2')
  }

  rowData(x@ProteinExperiment) <- value[[1]]
  rowData(x@PeptideExperiment) <- value[[2]]
  validObject(x)
  return(x)
})

#' @export
setGeneric('rowDataProt', function(x){
  standardGeneric('rowDataProt')
})

#' @export
setGeneric('rowDataProt<-', function(x, value){
  standardGeneric('rowDataProt<-')
})

#' @export
setMethod('rowDataProt', 'ProteinExperiment', function(x) {

  return(rowData(x))

})

#' @export
setMethod('rowDataProt<-', 'ProteinExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @export
setMethod('rowDataProt', 'ProteomicsExperiment', function(x) {

  return(rowData(x@ProteinExperiment))

})

#' @export
setMethod('rowDataProt<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @export
setGeneric('rowDataPept', function(x){
  standardGeneric('rowDataPept')
})

#' @export
setGeneric('rowDataPept<-', function(x, value){
  standardGeneric('rowDataPept<-')
})

#' @export
setMethod('rowDataPept', 'PeptideExperiment', function(x) {

  return(rowData(x))

})

#' @export
setMethod('rowDataPept<-', 'PeptideExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @export
setMethod('rowDataPept', 'ProteomicsExperiment', function(x) {

  return(rowData(x@PeptideExperiment))

})

#' @export
setMethod('rowDataPept<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
