#' @rdname classAccessors
#' @aliases assayNames
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assayNames
setMethod('assayNames', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assayNames<-
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assayNames<-
setMethod('assayNames<-',
          'ProteinExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assayNames
#' @usage NULL
#' @export
setMethod('assayNames', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assayNames<-
#' @usage NULL
#' @export
setMethod('assayNames<-',
          'PeptideExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assayNames
#' @usage NULL
#' @export
setMethod('assayNames', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assayNames(x@ProteinExperiment),
              peptide = assayNames(x@PeptideExperiment)))

})


#' @rdname classAccessors
#' @aliases assayNamesProt
#' @usage NULL
#' @export
setMethod('assayNamesProt', 'ProteinExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname classAccessors
#' @aliases `assayNamesProt<-`
#' @usage NULL
#' @export
setMethod('assayNamesProt<-', 'ProteinExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname classAccessors
#' @aliases assayNamesProt
#' @usage NULL
#' @export
setMethod('assayNamesProt', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@ProteinExperiment))

})

#' @rdname classAccessors
#' @aliases `assayNamesProt<-`
#' @usage NULL
#' @export
setMethod('assayNamesProt<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @aliases assayNamesPept
#' @usage NULL
#' @export
setMethod('assayNamesPept', 'PeptideExperiment', function(x) {

  return(assayNames(x))

})

#' @rdname classAccessors
#' @aliases `assayNamesPept<-`
#' @usage NULL
#' @export
setMethod('assayNamesPept<-', 'PeptideExperiment', function(x, value) {

  assayNames(x) <- value
  return(x)

})

#' @rdname classAccessors
#' @aliases assayNamesPept
#' @usage NULL
#' @export
setMethod('assayNamesPept', 'ProteomicsExperiment', function(x) {

  return(assayNames(x@PeptideExperiment))

})

#' @rdname classAccessors
#' @aliases `assayNamesPept<-`
#' @usage NULL
#' @export
setMethod('assayNamesPept<-', 'ProteomicsExperiment', function(x, value) {

  assayNames(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
