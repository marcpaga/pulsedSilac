#' @rdname classAccessors
#' @usage NULLs
#' @export
setMethod('metaoptions', 'ProteinExperiment', function(x){

  return(x@metaoptions)

})

#' @rdname classAccessors
#' @usage NULLs
#' @export
setMethod('metaoptions<-', 'ProteinExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @usage NULLs
#' @export
setMethod('metaoptions', 'PeptideExperiment', function(x){

  return(x@metaoptions)

})

#' @export
#' @rdname classAccessors
#' @usage NULL
setMethod('metaoptions<-', 'PeptideExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x)

})

#' @rdname classAccessors
#' @usage NULL
#' @export
setMethod('metaoptions', 'ProteomicsExperiment', function(x){

  return(x@metaoptions)

})

#' @rdname classAccessors
#' @usage NULL
setMethod('metaoptions<-', 'ProteomicsExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})

