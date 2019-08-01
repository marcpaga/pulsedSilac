#' @rdname classAccessors
#' @aliases rowData
#' @usage NULL
#' @export
setMethod('rowData', 'ProteinExperiment', function(x, use.names, ...){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases rowData<-
#' @usage NULL
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

#' @rdname classAccessors
#' @aliases rowData
#' @usage NULL
#' @export
setMethod('rowData', 'PeptideExperiment', function(x, use.names, ...){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases rowData<-
#' @usage NULL
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

#' @rdname classAccessors
#' @aliases rowData
#' @usage NULL
#' @export
setMethod('rowData', 'ProteomicsExperiment', function(x, use.names, ...){

  return(list(protein = rowData(x@ProteinExperiment),
              peptide = rowData(x@PeptideExperiment)))

})

#' @rdname classAccessors
#' @aliases rowData<-
#' @usage NULL
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


#' @rdname classAccessors
#' @aliases rowDataProt
#' @usage NULL
#' @export
setMethod('rowDataProt', 'ProteinExperiment', function(x) {

  return(rowData(x))

})

#' @rdname classAccessors
#' @aliases rowDataProt<-
#' @usage NULL
#' @export
setMethod('rowDataProt<-', 'ProteinExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname classAccessors
#' @aliases rowDataProt
#' @usage NULL
#' @export
setMethod('rowDataProt', 'ProteomicsExperiment', function(x) {

  return(rowData(x@ProteinExperiment))

})

#' @rdname classAccessors
#' @aliases rowDataProt<-
#' @usage NULL
#' @export
setMethod('rowDataProt<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname classAccessors
#' @aliases rowDataPept
#' @usage NULL
#' @export
setMethod('rowDataPept', 'PeptideExperiment', function(x) {

  return(rowData(x))

})

#' @rdname classAccessors
#' @aliases rowDataPept<-
#' @usage NULL
#' @export
setMethod('rowDataPept<-', 'PeptideExperiment', function(x, value) {

  rowData(x) <- value
  return(x)

})

#' @rdname classAccessors
#' @aliases rowDataPept
#' @usage NULL
#' @export
setMethod('rowDataPept', 'ProteomicsExperiment', function(x) {

  return(rowData(x@PeptideExperiment))

})

#' @rdname classAccessors
#' @aliases rowDataPept<-
#' @usage NULL
#' @export
setMethod('rowDataPept<-', 'ProteomicsExperiment', function(x, value) {

  rowData(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
