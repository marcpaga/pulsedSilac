#' @rdname SilacProteomicsExperiment-accessors
#' @aliases assayNames,SilacProteomicsExperiment-method
#' @export
setMethod('assayNames', 'SilacProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assayNames(x@SilacProteinExperiment),
              peptide = assayNames(x@SilacPeptideExperiment)))

})


#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesProt
#' @aliases assayNamesProt,SilacProteinExperiment-method
#' @export
setMethod('assayNamesProt', 'SilacProteinExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesProt<-
#' @aliases assayNamesProt<-,SilacProteinExperiment-method
#' @export
setMethod('assayNamesProt<-', 'SilacProteinExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesProt
#' @aliases assayNamesProt,SilacProteomicsExperiment-method
#' @export
setMethod('assayNamesProt', 'SilacProteomicsExperiment', function(x) {

  return(assayNames(x@SilacProteinExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesProt<-
#' @aliases assayNamesProt<-,SilacProteomicsExperiment-method
#' @export
setMethod('assayNamesProt<-', 'SilacProteomicsExperiment', function(x, value) {

  assayNames(x@SilacProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesPept
#' @aliases assayNamesPept,SilacPeptideExperiment-method
#' @export
setMethod('assayNamesPept', 'SilacPeptideExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesPept<-
#' @aliases assayNamesPept<-,SilacPeptideExperiment-method
#' @export
setMethod('assayNamesPept<-', 'SilacPeptideExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesPept
#' @aliases assayNamesPept,SilacProteomicsExperiment-method
#' @export
setMethod('assayNamesPept', 'SilacProteomicsExperiment', function(x) {

  return(assayNames(x@SilacPeptideExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name assayNamesPept<-
#' @aliases assayNamesPept<-,SilacProteomicsExperiment-method
#' @export
setMethod('assayNamesPept<-', 'SilacProteomicsExperiment', function(x, value) {

  assayNames(x@SilacPeptideExperiment) <- value
  validObject(x)
  return(x)

})
