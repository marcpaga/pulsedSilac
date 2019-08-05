#' @rdname ProteomicsExperiment-accessors
#' @aliases assayNames,ProteomicsExperiment-method
#' @export
setMethod('assayNames', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assayNames(x@ProteinExperiment),
              peptide = assayNames(x@PeptideExperiment)))

})


#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesProt
#' @aliases assayNamesProt,ProteinExperiment-method
#' @export
setMethod('assayNamesProt', 'ProteinExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesProt<-
#' @aliases assayNamesProt<-,ProteinExperiment-method
#' @export
setMethod('assayNamesProt<-', 'ProteinExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesProt
#' @aliases assayNamesProt,ProteomicsExperiment-method
#' @export
setMethod('assayNamesProt', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@ProteinExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesProt<-
#' @aliases assayNamesProt<-,ProteomicsExperiment-method
#' @export
setMethod('assayNamesProt<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesPept
#' @aliases assayNamesPept,PeptideExperiment-method
#' @export
setMethod('assayNamesPept', 'PeptideExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesPept<-
#' @aliases assayNamesPept<-,PeptideExperiment-method
#' @export
setMethod('assayNamesPept<-', 'PeptideExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesPept
#' @aliases assayNamesPept,ProteomicsExperiment-method
#' @export
setMethod('assayNamesPept', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@PeptideExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assayNamesPept<-
#' @aliases assayNamesPept<-,ProteomicsExperiment-method
#' @export
setMethod('assayNamesPept<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
