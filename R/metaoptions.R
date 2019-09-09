#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions,ProteinExperiment-method
#' @export
setMethod('metaoptions', 'ProteinExperiment', function(x){

  meta_opts <- metadata(x)[c('conditionCol', 'timeCol')]
  return(meta_opts)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,ProteinExperiment-method
#' @export
setMethod('metaoptions<-', 'ProteinExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  return(x)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions,PeptideExperiment-method
#' @export
setMethod('metaoptions', 'PeptideExperiment', function(x){

  meta_opts <- metadata(x)[c('conditionCol', 'timeCol', 'proteinCol')]
  return(meta_opts)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,PeptideExperiment-method
#' @export
setMethod('metaoptions<-', 'PeptideExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  return(x)

})

#' @rdname ProteomicsExperiment-accessors
#' @name metaoptions
#' @aliases metaoptions,ProteomicsExperiment-method
#' @export
setMethod('metaoptions', 'ProteomicsExperiment', function(x){

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'idColProt',
                         'idColPept',
                         'linkedSubset',
                         'subsetMode',
                         'proteinCol')
  meta_opts <- metadata(x)[metaoptions_names]
  return(meta_opts)

})

#' @rdname ProteomicsExperiment-accessors
#' @name metaoptions<-
#' @aliases metaoptions<-,ProteomicsExperiment-method
#' @export
setMethod('metaoptions<-', 'ProteomicsExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})

