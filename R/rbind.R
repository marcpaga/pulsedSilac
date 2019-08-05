#' @rdname ProteinPeptideExperiment-accessors
#' @aliases rbind,ProteinExperiment-method
#' @export
setMethod('rbind', 'ProteinExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  metaoptions(out) <- metaoptions(x)
  metadata(out) <- metadata(x)
  return(out)

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases rbind,PeptideExperiment-method
#' @export
setMethod('rbind', 'PeptideExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()

  metaoptions(out) <- metaoptions(x)
  metadata(out) <- metadata(x)
  return(out)

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases rbind,ProteomicsExperiment-method
#' @export
setMethod('rbind', 'ProteomicsExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  y <- unname(list(...))[[2]]

  new.prot <- rbind(x@ProteinExperiment, y@ProteinExperiment)
  new.pept <- rbind(x@PeptideExperiment, y@PeptideExperiment)

  new.linkerDf <- rbindLinkerDf(x@linkerDf, y@linkerDf)

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = new.prot,
            PeptideExperiment = new.pept,
            colData = x@colData,
            linkerDf = new.linkerDf,
            metadata = x@metadata,
            metaoptions = x@metaoptions)
  return(PE)

})
