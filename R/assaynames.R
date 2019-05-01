#' @export
#' @importFrom SummarizedExperiment assayNames
setMethod('assayNames', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
#' @importFrom SummarizedExperiment assayNames<-
setMethod('assayNames<-',
          'ProteinExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @export
setMethod('assayNames', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @export
setMethod('assayNames<-',
          'PeptideExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @export
setMethod('assayNames', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assayNames(x@ProteinExperiment),
              peptide = assayNames(x@PeptideExperiment)))

})


#' @export
setMethod('assayNamesProt', 'ProteinExperiment', function(x) {

  return(assayNames(x))

})

#' @export
setMethod('assayNamesProt<-', 'ProteinExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @export
setMethod('assayNamesProt', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@ProteinExperiment))

})

#' @export
setMethod('assayNamesProt<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @export
setMethod('assayNamesPept', 'PeptideExperiment', function(x) {

  return(assayNames(x))

})

#' @export
setMethod('assayNamesPept<-', 'PeptideExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @export
setMethod('assayNamesPept', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@PeptideExperiment))

})

#' @export
setMethod('assayNamesPept<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
