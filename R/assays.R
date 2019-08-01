#' @rdname classAccessors
#' @aliases assays
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assays
setMethod('assays', 'ProteinExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assays<-
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assays<-
setMethod('assays<-',
          'ProteinExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assays
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assays
setMethod('assays', 'PeptideExperiment', function(x, ..., withDimnames){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assays<-
#' @usage NULL
#' @export
#' @importFrom SummarizedExperiment assays<-
setMethod('assays<-',
          'PeptideExperiment', function(x, ..., withDimnames, value){

  return(callNextMethod())

})

#' @rdname classAccessors
#' @aliases assays
#' @usage NULL
#' @export
setMethod('assays', 'ProteomicsExperiment', function(x, ..., withDimnames){

  return(list(protein = assays(x@ProteinExperiment),
              peptide = assays(x@PeptideExperiment)))

})



#' @rdname classAccessors
#' @aliases assaysProt
#' @usage NULL
#' @export
setMethod('assaysProt', 'ProteinExperiment', function(x) {

  return(assays(x))

})

#' @rdname classAccessors
#' @aliases assaysProt<-
#' @usage NULL
#' @export
setMethod('assaysProt<-', 'ProteinExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @aliases assaysProt
#' @usage NULL
#' @export
setMethod('assaysProt', 'ProteomicsExperiment', function(x) {

  return(assays(x@ProteinExperiment))

})

#' @rdname classAccessors
#' @aliases assaysProt<-
#' @usage NULL
#' @export
setMethod('assaysProt<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@ProteinExperiment) <- value
  validObject(x)
  return(x)

})


#' @rdname classAccessors
#' @aliases assaysPept
#' @usage NULL
#' @export
setMethod('assaysPept', 'PeptideExperiment', function(x) {

  return(assays(x))

})

#' @rdname classAccessors
#' @aliases assaysPept<-
#' @usage NULL
#' @export
setMethod('assaysPept<-', 'PeptideExperiment', function(x, value) {

  assays(x) <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @aliases assaysPept
#' @usage NULL
#' @export
setMethod('assaysPept', 'ProteomicsExperiment', function(x) {

  return(assays(x@PeptideExperiment))

})

#' @rdname classAccessors
#' @aliases assaysPept<-
#' @usage NULL
#' @export
setMethod('assaysPept<-', 'ProteomicsExperiment', function(x, value) {

  assays(x@PeptideExperiment) <- value
  validObject(x)
  return(x)

})
