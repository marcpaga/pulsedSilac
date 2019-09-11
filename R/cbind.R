#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases cbind,SilacProteinExperiment-method
#' @export
setMethod('cbind', 'SilacProteinExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  out <- .removeDuplicatesMetaoptions(out)
  return(out)
})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases cbind,SilacPeptideExperiment-method
#' @export
setMethod('cbind', 'SilacPeptideExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  return(out)
})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases cbind,SilacProteomicsExperiment-method
#' @export
setMethod('cbind', 'SilacProteomicsExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  y <- unname(list(...))[[2]]

  new.prot <- cbind(x@SilacProteinExperiment, y@SilacProteinExperiment)
  new.pept <- cbind(x@SilacPeptideExperiment, y@SilacPeptideExperiment)

  PE <- new(Class = 'SilacProteomicsExperiment',
            SilacProteinExperiment = new.prot,
            SilacPeptideExperiment = new.pept,
            colData = colData(new.prot),
            linkerDf = x@linkerDf,
            metadata = x@metadata)
  return(PE)
})
