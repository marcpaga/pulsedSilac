#' @keywords internal
setMethod('show', 'ProteinExperiment', function(object){
  callNextMethod()
})

#' @keywords internal
setMethod('show', 'PeptideExperiment', function(object){
  callNextMethod()
})

#' Show method for ProteomicsExperiment
#'
#' @importMethodsFrom SummarizedExperiment show
#' @keywords internal
setMethod('show',
          'ProteomicsExperiment',
          function(object){

  selectSome <- S4Vectors:::selectSome
  scat <- function(fmt, vals=character(), exdent=2, startSpace = FALSE, ...)
  {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(S4Vectors:::selectSome(vals), collapse=" ")
    txt <- sprintf(fmt, length(vals), lbls)
    if (startSpace) {
      cat(paste0('  ',strwrap(txt, exdent=exdent, ...)), sep="\n")
    } else {
      cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }

  }

  ## class
  cat("class:", class(object), "\n")

  ## metadata
  expt <- names(metadata(object))
  if (is.null(expt)){
    expt <- character(length(metadata(object)))
  }
  scat("metadata(%d): %s\n", expt)

  ## protein part
  cat("protein:", "\n")
  cat("  |-- dim:", dim(object@ProteinExperiment), "\n")

  ## assaysProt
  nms <- assayNamesProt(object)
  if (is.null(nms)) {
    nms <- character(length(assaysProt(object, withDimnames=FALSE)))
  }
  scat("|-- assays(%d): %s\n", nms, startSpace = TRUE)

  ## rowDataProt
  dimnames <- rownames(rowDataProt(object))
  dlen <- sapply(dimnames, length)
  if (length(dlen) > 0) {
    scat("|-- rownames(%d): %s\n", dimnames[[1]], startSpace = TRUE)
  } else {
    scat("|-- rownames: NULL\n", startSpace = TRUE)
  }

  scat("|-- rowData names(%d): %s\n",
        names(rowDataProt(object)),  startSpace = TRUE)


  ## peptide part
  cat("peptide:", "\n")
  cat("  |-- dim:", dim(object@PeptideExperiment), "\n")

  nms <- assayNamesPept(object)
  if (is.null(nms)) {
    nms <- character(length(assaysPept(object, withDimnames=FALSE)))
  }
  scat("|-- assays(%d): %s\n", nms, startSpace = TRUE)

  ## rowDataPept
  dimnames <- rownames(rowDataPept(object))
  dlen <- sapply(dimnames, length)
  if (length(dlen) > 0) {
    scat("|-- rownames(%d): %s\n", dimnames[[1]], startSpace = TRUE)
  } else {
    scat("|-- rownames: NULL\n", startSpace = TRUE)
  }

  scat("|-- rowData names(%d): %s\n",
       names(rowDataPept(object)),  startSpace = TRUE)

  ## colData
  dimnames <- rownames(colData(object))
  if (length(dimnames)) scat("colnames(%d): %s\n", dimnames)
  else cat("colnames: NULL\n")
  scat("colData names(%d): %s\n", names(colData(object)), startSpace = FALSE)

})

testList[[3]]
