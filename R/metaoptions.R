#' @export
setMethod('metaoptions', 'ProteinExperiment', function(x){

  return(x@metaoptions)

})

#' @export
setMethod('metaoptions<-', 'ProteinExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x@metaoptions)

})

#' @export
setMethod('metaoptions', 'PeptideExperiment', function(x){

  return(x@metaoptions)

})

#' @export
setMethod('metaoptions<-', 'PeptideExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x@metaoptions)

})

#' @export
setMethod('metaoptions', 'ProteomicsExperiment', function(x){

  return(x@metaoptions)

})

setMethod('metaoptions<-', 'ProteomicsExperiment', function(x, value){

  x@metaoptions <- value
  validObject(x)
  return(x@metaoptions)

})

