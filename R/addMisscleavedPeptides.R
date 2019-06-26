#' @title Add miss-cleaved peptide expression levels for old isotope recycling
#' estimation
#'
#' @description To estimate how much of the "old" isotope is being used in
#' "new" proteins the expression level of miss-cleaved peptides that contain
#' a mix of isotopes (one old and one new) and miss-cleaved peptides that
#' contain only new isotopes can be used to estimate amino acid recycling.
#' To add this information  a \code{data.frame} with the following information
#' is required:
#' \itemize{
#'   \item A column with ids that can be mapped to the peptide rowData (not
#'   necessary for \code{ProteinExperiment} objects).
#'   \item A column indicating which isotope configuration the peptide has: two
#'   heavy isotopes, mix of light and heavy isotope, etc.
#'   \item A set of columns with peptide expression quantification. The number
#'   of columns should be the same, and in the same order, as the number of
#'   samples in the \code{ProteinExperiment}, \code{PeptideExperiment} or
#'   \code{ProteomicsExperiment}.
#' }
#' Non-detected measuremens should be \code{NA}.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment}.
#' @param newdata a \code{data.frame} containing the data described in the
#' description.
#' @param modCol \code{character} or \code{numeric} indicating which column of
#' newdata contains the peptide's isotope configuration.
#' @param dataCols \code{character} or \code{numeric} vector indicating which
#' columns of newdata contain the quantitative data. It should have the same
#' length as the number of samples in x (\code{ncol(x)}). Samples should be in
#' the same order.
#' @param idColPept \code{character} indicating which column contains the ids
#' that will be used in the merge with rowDataPept(x). The column in newdata
#' should have the same name as in rowDataPept(x).
#'
#' @return A \code{PeptideExperiment} or \code{ProteomicsExperiment} with new
#' assay entries. If x is a \code{ProteinExperiment} then a
#' \code{PeptideExperiment} is returned.
#'
#' @examples
#'
#' @export
setGeneric('addMisscleavedPeptides', function(x, ...){
  standardGeneric('addMisscleavedPeptides')
})


#' @export
setMethod('addMisscleavedPeptides',
          'ProteinExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {

  ## argument checks
  if (any(missing(x),
          missing(newdata),
          missing(modCol),
          missing(dataCols),
          missing(idColPept))) {
    stop('All arguments are mandatory')
  }

  if (!is.data.frame(newdata)) {
    stop('newdata must be a data.frame')
  }

  if (!is.numeric(modCol) & !is.character(modCol)) {
    stop('idCol must be numeric or character')
  }

  if (!is.numeric(dataCols) & !is.character(dataCols)) {
    stop('dataCols must be numeric or character')
  }

  if (!is.character(idColPept)) {
    stop('idColPept must be character')
  }

  ## get the amount of peptide configurations
  modifications <- unique(newdata[, modCol])

  ## remove peptides that are only in one modification
  toremove <- which(table(newdata[, idColPept]) < length(modifications))
  if (length(toremove) > 0) {
    newdata <- newdata[-which(newdata[, idColPept] %in% names(toremove)),]
  }

  ## for each configuration we will make an assay to put in assaysPep
  for (i in seq_along(modifications)) {

    ## for the assay name
    if (is.numeric(modCol)) {
      modColName <- colnames(newdata)[modCol]
    } else {
      modColName <- modCol
    }

    ## get the quantification data
    tempDf <- subset(newdata, get(modColName) == modifications[i])


    if (i == 1) {
      assaysList <- list()
      id_order <- as.vector(tempDf[, idColPept])
      rowData <- tempDf[, idColPept, drop = FALSE]
    } else {
      tempDf <- tempDf[match(tempDf[, idColPept], id_order),]
    }

    assayMat <- as.matrix(tempDf[, dataCols])
    rownames(assayMat) <- NULL
    colnames(assayMat) <- NULL
    assaysList[[i]] <- assayMat
  }

  names(assaysList) <- modifications


  pe <- PeptideExperiment(assays = assaysList,
                    rowData = rowData,
                    colData = colData(x),
                    conditionCol = metaoptions(x)[['conditionCol']],
                    timeCol = metaoptions(x)[['timeCol']],
                    replicateIntCol = metaoptions(x)[['replicateIntCol']],
                    replicateTimeCol = metaoptions(x)[['replicateTimeCol']])

  return(pe)
})

#' @export
setMethod('addMisscleavedPeptides',
          'PeptideExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {

  ## argument checks
  if (any(missing(x),
          missing(newdata),
          missing(modCol),
          missing(dataCols),
          missing(idColPept))) {
    stop('All arguments are mandatory')
  }

  if (!is.data.frame(newdata)) {
    stop('newdata must be a data.frame')
  }

  if (!is.numeric(modCol) & !is.character(modCol)) {
    stop('idCol must be numeric or character')
  }

  if (!is.numeric(dataCols) & !is.character(dataCols)) {
    stop('dataCols must be numeric or character')
  }

  if (!is.character(idColPept)) {
    stop('idColPept must be character')
  }

  ## get the amount of peptide configurations
  modifications <- unique(newdata[, modCol])

  ## for each configuration we will make an assay to put in assaysPep
  for (i in seq_along(modifications)) {

    ## for the assay name
    if (is.numeric(modCol)) {
      modColName <- colnames(newdata)[modCol]
    } else {
      modColName <- modCol
    }

    ## get the quantification data
    tempDf <- subset(newdata, get(modColName) == modifications[i])
    dataMat <- as.matrix(tempDf[, dataCols])
    colnames(dataMat) <- NULL
    rownames(dataMat) <- NULL

    ## make the assay-like matrix, we use try in case no assayPep are present
    ## yets
    dimPep <- try(dim(assays(x)[[1]]))
    if (is(dimPep, 'try-error')) {
      dimPep <- c(nrow(rowData(x)), nrow(colData(x)))
    }
    tempMat <- matrix(NA, nrow = dimPep[1], ncol = dimPep[2])

    ## fill the new assay with the data using the match
    newrows <- match(tempDf[, idColPept], rowData(x)[,idColPept])
    oldrows <- seq_len(nrow(dataMat))
    nomatch <- which(is.na(newrows))

    ## remove peptides that cant be matched
    oldrows <- oldrows[-nomatch]
    newrows <- newrows[-nomatch]

    if (length(oldrows) == 0) {
      warning('None of the given ids were matched')
    }

    tempMat[newrows, ] <- dataMat[oldrows, ]
    assays(x)[[as.character(modifications[i])]] <- tempMat

  }

  return(x)

})

#' @export
setMethod('addMisscleavedPeptides',
          'ProteomicsExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {


  new.pept <- addMisscleavedPeptides(x = x@PeptideExperiment,
                                    newdata = newdata,
                                    modCol = modCol,
                                    dataCols = dataCols,
                                    idColPept = idColPept)

  x@PeptideExperiment <- new.pept
  return(x)
})

