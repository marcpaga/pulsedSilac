#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions,ProteinExperiment-method
#' @export
setMethod('metaoptions', 'ProteinExperiment', function(x){

  return(x@metaoptions)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,ProteinExperiment-method
#' @export
setMethod('metaoptions<-', 'ProteinExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions,PeptideExperiment-method
#' @export
setMethod('metaoptions', 'PeptideExperiment', function(x){

  return(x@metaoptions)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,PeptideExperiment-method
#' @export
setMethod('metaoptions<-', 'PeptideExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name metaoptions
#' @aliases metaoptions,ProteomicsExperiment-method
#' @export
setMethod('metaoptions', 'ProteomicsExperiment', function(x){

  return(x@metaoptions)

})

#' @rdname ProteomicsExperiment-accessors
#' @name metaoptions<-
#' @aliases metaoptions<-,ProteomicsExperiment-method
#' @export
setMethod('metaoptions<-', 'ProteomicsExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})

