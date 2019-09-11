#' @rdname SilacProteinPeptideExperiment-accessors
#' @name SilacProteinPeptideExperiment-accessors
#' @aliases assays<-,SilacProteinExperiment,ANY-method
#' @export
setMethod('assays<-', 'SilacProteinExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @name SilacProteinPeptideExperiment-accessors
#' @aliases assays<-,SilacPeptideExperiment,ANY-method
#' @export
setMethod('assays<-', 'SilacPeptideExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assays
#' @aliases assays,SilacProteomicsExperiment-method
#' @export
setMethod('assays', 'SilacProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assays(x@SilacProteinExperiment),
              peptide = assays(x@SilacPeptideExperiment)))

})


#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysProt
#' @aliases assaysProt,SilacProteinExperiment-method
#' @export
setMethod('assaysProt', 'SilacProteinExperiment', function(x) {

  return(assays(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysProt<-
#' @aliases assaysProt<-,SilacProteinExperiment-method
#' @export
setMethod('assaysProt<-', 'SilacProteinExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysProt
#' @aliases assaysProt,SilacProteomicsExperiment-method
#' @export
setMethod('assaysProt', 'SilacProteomicsExperiment', function(x) {

  return(assays(x@SilacProteinExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysProt<-
#' @aliases assaysProt<-,SilacProteomicsExperiment-method
#' @export
setMethod('assaysProt<-', 'SilacProteomicsExperiment', function(x, value) {

  assays(x@SilacProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysPept
#' @aliases assaysPept,SilacPeptideExperiment-method
#' @export
setMethod('assaysPept', 'SilacPeptideExperiment', function(x) {

  return(assays(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysPept<-
#' @aliases assaysPept<-,SilacPeptideExperiment-method
#' @export
setMethod('assaysPept<-', 'SilacPeptideExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysPept
#' @aliases assaysPept,SilacProteomicsExperiment-method
#' @export
setMethod('assaysPept', 'SilacProteomicsExperiment', function(x) {

  return(assays(x@SilacPeptideExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assaysPept<-
#' @aliases assaysPept<-,SilacProteomicsExperiment-method
#' @export
setMethod('assaysPept<-', 'SilacProteomicsExperiment', function(x, value) {

  assays(x@SilacPeptideExperiment) <- value
  validObject(x)
  return(x)

})
