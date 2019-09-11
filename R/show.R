#' @keywords internal
setMethod('show', 'SilacProteinExperiment', function(object){
  callNextMethod()
})

#' @keywords internal
setMethod('show', 'SilacPeptideExperiment', function(object){
  callNextMethod()
})

#' @importMethodsFrom SummarizedExperiment show
#' @keywords internal
setMethod('show',
          'SilacProteomicsExperiment',
          function(object){

  ## this is the unexported function from the S4Vectors:::selectSome
  selectSome <- function (obj, maxToShow = 5, ellipsis = "...",
                          ellipsisPos = c("middle", "end", "start"),
                          quote = FALSE) {
    if (is.character(obj) && quote)
      obj <- sQuote(obj)
    ellipsisPos <- match.arg(ellipsisPos)
    len <- length(obj)
    if (maxToShow < 3)
      maxToShow <- 3
    if (len > maxToShow) {
      maxToShow <- maxToShow - 1
      if (ellipsisPos == "end") {
        c(head(obj, maxToShow), ellipsis)
      }
      else if (ellipsisPos == "start") {
        c(ellipsis, tail(obj, maxToShow))
      }
      else {
        bot <- ceiling(maxToShow/2)
        top <- len - (maxToShow - bot - 1)
        nms <- obj[c(seq_len(bot), top:len)]
        c(as.character(nms[seq_len(bot)]), ellipsis, as.character(nms[-c(seq_len(bot))]))
      }
    }
    else {
      obj
    }
  }

  ## copied from the SummarizedExperiment show method
  scat <- function(fmt, vals=character(), exdent=2, startSpace = FALSE, ...)
  {
    vals <- ifelse(nzchar(vals), vals, "''")
    lbls <- paste(selectSome(vals), collapse=" ")
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
  cat("  |-- dim:", dim(object@SilacProteinExperiment), "\n")

  ## assaysProt
  nms <- assayNamesProt(object)
  if (is.null(nms)) {
    nms <- character(length(assaysProt(object)))
  }
  scat("|-- assays(%d): %s\n", nms, startSpace = TRUE)

  ## rowDataProt
  dimnames <- rownames(rowDataProt(object))
  dlen <- vapply(dimnames, FUN.VALUE = 1, FUN = length)
  if (length(dlen) > 0) {
    scat("|-- rownames(%d): %s\n", dimnames, startSpace = TRUE)
  } else {
    scat("|-- rownames: NULL\n", startSpace = TRUE)
  }

  scat("|-- rowData names(%d): %s\n",
        names(rowDataProt(object)),  startSpace = TRUE)


  ## peptide part
  cat("peptide:", "\n")
  cat("  |-- dim:", dim(object@SilacPeptideExperiment), "\n")

  nms <- assayNamesPept(object)
  if (is.null(nms)) {
    nms <- character(length(assaysPept(object)))
  }
  scat("|-- assays(%d): %s\n", nms, startSpace = TRUE)

  ## rowDataPept
  dimnames <- rownames(rowDataPept(object))
  dlen <- vapply(dimnames, FUN.VALUE = 1, FUN = length)
  if (length(dlen) > 0) {
    scat("|-- rownames(%d): %s\n", dimnames, startSpace = TRUE)
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

