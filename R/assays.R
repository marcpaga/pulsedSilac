#' @export
setMethod('assays', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
setMethod('assays<-',
          'ProteinExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @export
setMethod('assays', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
setMethod('assays<-',
          'PeptideExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @export
setMethod('assays', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assays(x@ProteinExperiment),
              peptide = assays(x@PeptideExperiment)))

})

#' @export
setMethod('assays<-',
          'ProteomicsExperiment', function(x, ..., withDimnames, value){

  if (!is.list(value) | length(value) != 2) {
    stop('Value must be a list of length 2')
  }

  assays(x@ProteinExperiment) <- value[[1]]
  assays(x@PeptideExperiment) <- value[[2]]
  validObject(x)
  return(x)
})

#' @export
setGeneric('assaysProt', function(x){
  standardGeneric('assaysProt')
})

#' @export
setGeneric('assaysProt<-', function(x, value){
  standardGeneric('assaysProt<-')
})

#' @export
setMethod('assaysProt', 'ProteinExperiment', function(x) {

  return(assays(x))

})

#' @export
setMethod('assaysProt<-', 'ProteinExperiment', function(x, value) {

  assays(x) <- value
  return(x)

})

#' @export
setMethod('assaysProt', 'ProteomicsExperiment', function(x) {

  return(assays(x@ProteinExperiment))

})

#' @export
setMethod('assaysProt<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @export
setGeneric('assaysPept', function(x){
  standardGeneric('assaysPept')
})

#' @export
setGeneric('assaysPept<-', function(x, value){
  standardGeneric('assaysPept<-')
})

#' @export
setMethod('assaysPept', 'PeptideExperiment', function(x) {

  return(assays(x))

})

#' @export
setMethod('assaysPept<-', 'PeptideExperiment', function(x, value) {

  assays(x) <- value
  return(x)

})

#' @export
setMethod('assaysPept', 'ProteomicsExperiment', function(x) {

  return(assays(x@PeptideExperiment))

})

#' @export
setMethod('assaysPept<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
