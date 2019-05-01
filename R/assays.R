#' @export
#' @importFrom SummarizedExperiment assays
setMethod('assays', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
#' @importFrom SummarizedExperiment assays<-
setMethod('assays<-',
          'ProteinExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @export
#' @importFrom SummarizedExperiment assays
setMethod('assays', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
#' @importFrom SummarizedExperiment assays<-
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
