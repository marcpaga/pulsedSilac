#' @rdname filterByMissingTimepoints
#' @name filterByMissingTimepoints
#' @title Filter proteins/peptides by the amount of measurements overtime
#'
#' @description Searches for proteins/peptides that are not found in all
#' timepoints. This can be done for each condition independently
#' (strict = FALSE) or shared across conditions (strict = TRUE).
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param assayName A character indicating which assay will be used to count
#' the number of missed timepoints.
#' @param maxMissing A numeric indicating how many timepoints are allowed to be
#' missed.
#' @param strict Logical: if TRUE, then proteins have to meet the maxMissing
#' criteria in all conditions and time replicates to pass; if FALSE then
#' proteins only have to meet the maxMissing criteria in one condition or time
#' replicate to pass.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param replicateTimeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different time replicates.
#' @param returnVector Logical: if TRUE then a vector with the positions to be
#' subset is returned.
#' @param ... Unused.
#'
#' @return A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object or a logical vector with the rows that
#' pass the minimum number of desired timepoints.
#'
#' @examples
#' filterByMissingTimepoints(wormsPE,
#'                           assayName = 'ratio',
#'                           maxMissing = 2,
#'                           strict = FALSE)
#' @export
setGeneric('filterByMissingTimepoints', function(x, ...){
  standardGeneric('filterByMissingTimepoints')
})

#' @rdname filterByMissingTimepoints
#' @export
setMethod('filterByMissingTimepoints',
          'ProteinExperiment',
          function(x,
                   assayName,
                   maxMissing = 0,
                   strict = TRUE,
                   conditionCol,
                   replicateTimeCol,
                   returnVector = FALSE) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## assay used for the plotting
  mat <- assays(x)[[assayName]]

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  ## first it tries to get both condition and time replicates, if it fails then
  ## only condition and if it fails then all the samples are grouped together
  loopCols <- tryCatch(
    {
      loopCols <- experimentLoopWrapper(x, 'cond.timerep')
    },
    error = function(c){
      tryCatch(
        {
          loopCols <- experimentLoopWrapper(x, 'cond')
        },
        error = function(c){
          loopCols <- list(seq_len(ncol(x)))
        }
      )
    }
  )

  for (i in seq_along(loopCols)) {
    if (i == 1) {
      logMatrix <- matrix(data = FALSE,
                          ncol = length(loopCols),
                          nrow = nrow(x))
    }

    nacounts <- apply(mat[, loopCols[[i]]], 1, function(x) sum(is.na(x)))
    ## TRUE if there are less or equal missing as maxMissing
    logMatrix[, i] <- (nacounts <= maxMissing)

  }

  if (strict) {
    subsetVec <- apply(logMatrix, 1, all)
  } else {
    subsetVec <- apply(logMatrix, 1, any)
  }

  if (returnVector) {
    return(subsetVec)
  }

  new.x <- x[which(subsetVec),]
  validObject(new.x)
  return(new.x)

})


#' @rdname filterByMissingTimepoints
#' @export
setMethod('filterByMissingTimepoints',
          'PeptideExperiment',
          function(x,
                   assayName,
                   maxMissing = 0,
                   strict = TRUE,
                   conditionCol,
                   replicateTimeCol,
                   returnVector = FALSE) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## assay used for the plotting
  mat <- assays(x)[[assayName]]

  ## if the metaoptions are not given, try to get them from the slot
  ## if there are also not present, then the plot does not take into account
  ## conditions.
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }

  ## use trycatch since giveMetaoption raises and error if it does not find it,
  ## but for plotting metaoptions are not strictly necessary
  ## first it tries to get both condition and time replicates, if it fails then
  ## only condition and if it fails then all the samples are grouped together
  loopCols <- tryCatch(
    {
      loopCols <- experimentLoopWrapper(x, 'cond.timerep')
    },
    error = function(c){
      tryCatch(
        {
          loopCols <- experimentLoopWrapper(x, 'cond')
        },
        error = function(c){
          loopCols <- list(seq_len(ncol(x)))
        }
      )
    }
  )

  for (i in seq_along(loopCols)) {
    if (i == 1) {
      logMatrix <- matrix(data = FALSE,
                          ncol = length(loopCols),
                          nrow = nrow(x))
    }

    nacounts <- apply(mat[, loopCols[[i]]], 1, function(x) sum(is.na(x)))
    logMatrix[, i] <- (nacounts <= maxMissing)

  }

  if (strict) {
    subsetVec <- apply(logMatrix, 1, all)
  } else {
    subsetVec <- apply(logMatrix, 1, any)
  }

  if (returnVector) {
    return(subsetVec)
  }

  new.x <- x[which(subsetVec),]
  validObject(new.x)
  return(new.x)

})

#' @rdname filterByMissingTimepoints
#' @export
setMethod('filterByMissingTimepoints',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   maxMissing = 0,
                   strict = TRUE,
                   conditionCol,
                   replicateTimeCol,
                   returnVector = FALSE) {

  if (giveMetaoption(x, 'subsetMode') == 'peptide') {

    subsetVec <- filterByMissingTimepoints(x@PeptideExperiment,
                                           assayName = assayName,
                                           maxMissing = maxMissing,
                                           strict = strict,
                                           conditionCol,
                                           replicateTimeCol,
                                           returnVector = TRUE)

  } else {

    subsetVec <- filterByMissingTimepoints(x@ProteinExperiment,
                                           assayName = assayName,
                                           maxMissing = maxMissing,
                                           strict = strict,
                                           conditionCol,
                                           replicateTimeCol,
                                           returnVector = TRUE)
  }

  if (returnVector) {
    return(subsetVec)
  }

  new.x <- x[which(subsetVec),]
  validObject(new.x)
  return(new.x)

})
