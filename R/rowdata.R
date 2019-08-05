#' @rdname ProteinPeptideExperiment-accessors
#' @aliases rowData<-,ProteinExperiment-method
#' @importFrom S4Vectors DataFrame
#' @export
setMethod('rowData<-', 'ProteinExperiment', function(x, ..., value){

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  if (is(value, 'DataFrame')) {
    x@elementMetadata <- value
  }

  validObject(x)
  return(x)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases rowData<-,PeptideExperiment-method
#' @importFrom S4Vectors DataFrame
#' @export
setMethod('rowData<-', 'PeptideExperiment', function(x, ..., value){

  if (is(value, 'data.frame')) {
    value <- DataFrame(value)
  }

  if (is(value, 'DataFrame')) {
    x@elementMetadata <- value
  }

  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases rowData,ProteomicsExperiment-method
#' @export
setMethod('rowData', 'ProteomicsExperiment', function(x, use.names, ...){

  return(list(protein = rowData(x@ProteinExperiment),
              peptide = rowData(x@PeptideExperiment)))

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases rowData<-,ProteomicsExperiment-method
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


#' @rdname ProteomicsExperiment-accessors
#' @name rowDataProt
#' @aliases rowDataProt,ProteinExperiment-method
#' @export
setMethod('rowDataProt', 'ProteinExperiment', function(x) {

  return(rowData(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataProt<-
#' @aliases rowDataProt<-,ProteinExperiment-method
#' @export
setMethod('rowDataProt<-', 'ProteinExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataProt
#' @aliases rowDataProt,ProteomicsExperiment-method
#' @export
setMethod('rowDataProt', 'ProteomicsExperiment', function(x) {

  return(rowData(x@ProteinExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataProt<-
#' @aliases rowDataProt<-,ProteomicsExperiment-method
#' @export
setMethod('rowDataProt<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname ProteomicsExperiment-accessors
#' @name rowDataPept
#' @aliases rowDataPept,PeptideExperiment-method
#' @export
setMethod('rowDataPept', 'PeptideExperiment', function(x) {

  return(rowData(x))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataPept<-
#' @aliases rowDataPept<-,PeptideExperiment-method
#' @export
setMethod('rowDataPept<-', 'PeptideExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataPept
#' @aliases rowDataPept,ProteomicsExperiment-method
#' @export
setMethod('rowDataPept', 'ProteomicsExperiment', function(x) {

  return(rowData(x@PeptideExperiment))

})

#' @rdname ProteomicsExperiment-accessors
#' @name rowDataPept<-
#' @aliases rowDataPept<-,ProteomicsExperiment-method
#' @usage NULL
#' @export
setMethod('rowDataPept<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
