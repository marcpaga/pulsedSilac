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
#'   necessary for \code{SilacProteinExperiment} objects).
#'   \item A column indicating which isotope configuration the peptide has: two
#'   heavy isotopes, mix of light and heavy isotope, etc.
#'   \item A set of columns with peptide expression quantification. The number
#'   of columns should be the same, and in the same order, as the number of
#'   samples in the \code{SilacProteinExperiment}, \code{SilacPeptideExperiment}
#'   or \code{SilacProteomicsExperiment}.
#' }
#' Non-detected measurements should be \code{NA}.
#'
#' @param x A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or
#' \code{SilacProteomicsExperiment}.
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
#' @param ... Unused.
#'
#' @return A \code{SilacPeptideExperiment} or \code{SilacProteomicsExperiment}
#' with new assay entries. If x is a \code{SilacProteinExperiment} then a
#' \code{SilacPeptideExperiment} is returned.
#'
#' @examples
#' data('wormsPE')
#' data('recycleLightLysine')
#' protPE <- ProtExp(wormsPE)
#' missPE <- addMisscleavedPeptides(x = protPE,
#'                                  newdata = recycleLightLysine,
#'                                  idColPept = 'Sequence',
#'                                  modCol = 'Modifications',
#'                                  dataCols = c(18:31))
#'
#' names(assays(missPE))[1:2] <- c('int_lys8lys8', 'int_lys8lys0')
#' missPE <- calculateOldIsotopePool(x = missPE, 'int_lys8lys8', 'int_lys8lys0')
#'
#' plotDistributionAssay(missPE, assayName = 'oldIsotopePool')
#'
#' @export
setGeneric('addMisscleavedPeptides', function(x, ...){
  standardGeneric('addMisscleavedPeptides')
})


#' @rdname addMisscleavedPeptides
#' @export
setMethod('addMisscleavedPeptides',
          'SilacProteinExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {

  ## argument checks -----------------------------------------------------------
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

  ## data processing -----------------------------------------------------------
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


  pe <- SilacPeptideExperiment(assays = assaysList,
                    rowData = rowData,
                    colData = colData(x),
                    conditionCol = metaoptions(x)[['conditionCol']],
                    timeCol = metaoptions(x)[['timeCol']])

  return(pe)
})

#' @rdname addMisscleavedPeptides
#' @export
setMethod('addMisscleavedPeptides',
          'SilacPeptideExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {

  ## argument checker ----------------------------------------------------------
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

  ## data processing -----------------------------------------------------------
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
    newrows <- match(tempDf[, idColPept], rowData(x)[, idColPept])
    oldrows <- seq_len(nrow(dataMat))
    nomatch <- which(is.na(newrows))

    ## remove peptides that cant be matched
    if (length(nomatch) > 0) {
      oldrows <- oldrows[-nomatch]
      newrows <- newrows[-nomatch]
    }

    if (length(oldrows) == 0) {
      warning('None of the given ids were matched')
    }

    tempMat[newrows, ] <- dataMat[oldrows, ]
    assays(x)[[as.character(modifications[i])]] <- tempMat

  }

  return(x)

})

#' @rdname addMisscleavedPeptides
#' @export
setMethod('addMisscleavedPeptides',
          'SilacProteomicsExperiment',
          function(x,
                   newdata,
                   modCol,
                   dataCols,
                   idColPept) {


  new.pept <- addMisscleavedPeptides(x = x@SilacPeptideExperiment,
                                    newdata = newdata,
                                    modCol = modCol,
                                    dataCols = dataCols,
                                    idColPept = idColPept)

  x@SilacPeptideExperiment <- new.pept
  return(x)
})

