#' @export
setMethod('cbind', 'ProteinExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()
  metaoptions(out) <- metaoptions(x)
  return(out)
})

#' @export
setMethod('cbind', 'PeptideExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  out <- callNextMethod()
  metaoptions(out) <- metaoptions(x)
  return(out)
})

#' @export
setMethod('cbind', 'ProteomicsExperiment', function(..., deparse.level = 1) {

  x <- unname(list(...))[[1]]
  y <- unname(list(...))[[2]]

  new.prot <- cbind(x@ProteinExperiment, y@ProteinExperiment)
  new.pept <- cbind(x@PeptideExperiment, y@PeptideExperiment)

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = new.prot,
            PeptideExperiment = new.pept,
            colData = colData(new.prot),
            linkerDf = x@linkerDf,
            metadata = x@metadata,
            metaoptions = x@metaoptions)
  return(PE)
})
