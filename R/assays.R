#' @rdname ProteinPeptideExperiment-accessors
#' @name ProteinPeptideExperiment-accessors
#' @aliases assays<-,ProteinExperiment,ANY-method
#' @export
setMethod('assays<-', 'ProteinExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @name ProteinPeptideExperiment-accessors
#' @aliases assays<-,PeptideExperiment,ANY-method
#' @export
setMethod('assays<-', 'PeptideExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assays
#' @aliases assays,ProteomicsExperiment-method
#' @export
setMethod('assays', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assays(x@ProteinExperiment),
              peptide = assays(x@PeptideExperiment)))

})


#' @rdname ProteomicsExperiment-accessors
#' @name assaysProt
#' @aliases assaysProt,ProteinExperiment-method
#' @export
setMethod('assaysProt', 'ProteinExperiment', function(x) {

  return(assays(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysProt<-
#' @aliases assaysProt<-,ProteinExperiment-method
#' @export
setMethod('assaysProt<-', 'ProteinExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysProt
#' @aliases assaysProt,ProteomicsExperiment-method
#' @export
setMethod('assaysProt', 'ProteomicsExperiment', function(x) {

  return(assays(x@ProteinExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysProt<-
#' @aliases assaysProt<-,ProteomicsExperiment-method
#' @export
setMethod('assaysProt<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname ProteomicsExperiment-accessors
#' @name assaysPept
#' @aliases assaysPept,PeptideExperiment-method
#' @export
setMethod('assaysPept', 'PeptideExperiment', function(x) {

  return(assays(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysPept<-
#' @aliases assaysPept<-,PeptideExperiment-method
#' @export
setMethod('assaysPept<-', 'PeptideExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysPept
#' @aliases assaysPept,ProteomicsExperiment-method
#' @export
setMethod('assaysPept', 'ProteomicsExperiment', function(x) {

  return(assays(x@PeptideExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name assaysPept<-
#' @aliases assaysPept<-,ProteomicsExperiment-method
#' @export
setMethod('assaysPept<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
