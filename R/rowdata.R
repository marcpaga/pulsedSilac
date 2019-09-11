#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases rowData<-,SilacProteinExperiment-method
#' @importFrom S4Vectors DataFrame
#' @export
setMethod('rowData<-', 'SilacProteinExperiment', function(x, ..., value){

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  if (is(value, 'DataFrame')) {
    x@elementMetadata <- value
  }

  validObject(x)
  return(x)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases rowData<-,SilacPeptideExperiment-method
#' @importFrom S4Vectors DataFrame
#' @export
setMethod('rowData<-', 'SilacPeptideExperiment', function(x, ..., value){

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  if (is(value, 'DataFrame')) {
    x@elementMetadata <- value
  }

  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases rowData,SilacProteomicsExperiment-method
#' @export
setMethod('rowData', 'SilacProteomicsExperiment', function(x, use.names, ...){

  return(list(protein = rowData(x@SilacProteinExperiment),
              peptide = rowData(x@SilacPeptideExperiment)))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases rowData<-,SilacProteomicsExperiment-method
#' @export
setMethod('rowData<-', 'SilacProteomicsExperiment', function(x, ..., value){

  if (!is.list(value) | length(value) != 2) {
    stop('Value must be a list of length 2')
  }

  rowData(x@SilacProteinExperiment) <- value[[1]]
  rowData(x@SilacPeptideExperiment) <- value[[2]]
  validObject(x)
  return(x)
})


#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataProt
#' @aliases rowDataProt,SilacProteinExperiment-method
#' @export
setMethod('rowDataProt', 'SilacProteinExperiment', function(x) {

  return(rowData(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataProt<-
#' @aliases rowDataProt<-,SilacProteinExperiment-method
#' @export
setMethod('rowDataProt<-', 'SilacProteinExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataProt
#' @aliases rowDataProt,SilacProteomicsExperiment-method
#' @export
setMethod('rowDataProt', 'SilacProteomicsExperiment', function(x) {

  return(rowData(x@SilacProteinExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataProt<-
#' @aliases rowDataProt<-,SilacProteomicsExperiment-method
#' @export
setMethod('rowDataProt<-', 'SilacProteomicsExperiment', function(x, value) {

  rowData(x@SilacProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataPept
#' @aliases rowDataPept,SilacPeptideExperiment-method
#' @export
setMethod('rowDataPept', 'SilacPeptideExperiment', function(x) {

  return(rowData(x))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataPept<-
#' @aliases rowDataPept<-,SilacPeptideExperiment-method
#' @export
setMethod('rowDataPept<-', 'SilacPeptideExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataPept
#' @aliases rowDataPept,SilacProteomicsExperiment-method
#' @export
setMethod('rowDataPept', 'SilacProteomicsExperiment', function(x) {

  return(rowData(x@SilacPeptideExperiment))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name rowDataPept<-
#' @aliases rowDataPept<-,SilacProteomicsExperiment-method
#' @usage NULL
#' @export
setMethod('rowDataPept<-', 'SilacProteomicsExperiment', function(x, value) {

  rowData(x@SilacPeptideExperiment) <- value
  validObject(x)
  return(x)

})
