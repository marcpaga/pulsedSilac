#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases rbind,SilacProteinExperiment-method
#' @export
setMethod('rbind', 'SilacProteinExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  metaoptions(out) <- metaoptions(x)
  metadata(out) <- metadata(x)
  return(out)

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases rbind,SilacPeptideExperiment-method
#' @export
setMethod('rbind', 'SilacPeptideExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  metaoptions(out) <- metaoptions(x)
  metadata(out) <- metadata(x)
  return(out)

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases rbind,SilacProteomicsExperiment-method
#' @export
setMethod('rbind', 'SilacProteomicsExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  y <- unname(list(...))[[2]]

  new.prot <- rbind(x@SilacProteinExperiment, y@SilacProteinExperiment)
  new.pept <- rbind(x@SilacPeptideExperiment, y@SilacPeptideExperiment)

  new.linkerDf <- rbindLinkerDf(x@linkerDf, y@linkerDf)

  PE <- new(Class = 'SilacProteomicsExperiment',
            SilacProteinExperiment = new.prot,
            SilacPeptideExperiment = new.pept,
            colData = x@colData,
            linkerDf = new.linkerDf,
            metadata = x@metadata)
  return(PE)

})
