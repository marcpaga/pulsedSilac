#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases metaoptions,SilacProteinExperiment-method
#' @export
setMethod('metaoptions', 'SilacProteinExperiment', function(x){

  meta_opts <- metadata(x)[c('conditionCol', 'timeCol')]
  return(meta_opts)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,SilacProteinExperiment-method
#' @export
setMethod('metaoptions<-', 'SilacProteinExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases metaoptions,SilacPeptideExperiment-method
#' @export
setMethod('metaoptions', 'SilacPeptideExperiment', function(x){

  meta_opts <- metadata(x)[c('conditionCol', 'timeCol', 'proteinCol')]
  return(meta_opts)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases metaoptions<-,SilacPeptideExperiment-method
#' @export
setMethod('metaoptions<-', 'SilacPeptideExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  return(x)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name metaoptions
#' @aliases metaoptions,SilacProteomicsExperiment-method
#' @export
setMethod('metaoptions', 'SilacProteomicsExperiment', function(x){

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

#' @rdname SilacProteomicsExperiment-accessors
#' @name metaoptions<-
#' @aliases metaoptions<-,SilacProteomicsExperiment-method
#' @export
setMethod('metaoptions<-', 'SilacProteomicsExperiment', function(x, value){

  x@metadata <- value
  validObject(x)
  x <- synchronizeMetaoptions(x)
  return(x)

})

