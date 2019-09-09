#' @rdname scatterCompareAssays
#' @name scatterCompareAssays
#' @title Scatter plot of two conditions for each timepoint of an assay.
#'
#' @description Scatter plot of two conditions/replicates for a selected assay.
#' Timepoints are separated using \code{facet_wrap}.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param conditions A \code{character} of length 2 indicating which 2
#' conditions should be compared.
#' @param assayName Name of the assay to use in the plot.
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#' @param mode A \code{character} indicating which level of data to use,
#' either "protein" or "peptide". Only relevant for ProteomicsExperiment
#' inputs.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param timeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different timepoints.
#' @param ... Unused.
#'
#' @return A \code{ggplot} object or the \code{data.frame} that would be used
#' instead in the plot.
#'
#' @examples
#' data('wormsPE')
#' scatterCompareAssays(x = wormsPE,
#'                      conditions = c('OW40', 'OW450'),
#'                      assayName = 'ratio',
#'                      mode = 'protein')
#'
#' @import ggplot2
#' @export
setGeneric('scatterCompareAssays', function(x, ...){
  standardGeneric('scatterCompareAssays')
})

#' @rdname scatterCompareAssays
#' @export
setMethod('scatterCompareAssays', 'ProteinExperiment',
          function(x,
                   conditions,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {


  ## argument checker ----------------------------------------------------------
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }
  if (!missing(conditionCol)) {
    metadata(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metadata(x)[['timeCol']] <- timeCol
  }
  if (length(conditions) != 2) {
    stop('conditions must be a character vector of lenght 2')
  }

  ## data processing -----------------------------------------------------------
  ## extract the matrix
  mat <- assays(x)[[assayName]]

  loopCols <- .loopWrapper(x, 'conditionCol')
  if (length(loopCols) < 2) {
    stop('There is only one condition, comparisons cannot be made')
  }

  ## reduce loopCols to the two conditions
  availableConditions <- names(loopCols)
  if (!all(conditions %in% availableConditions)) {
    txt <- c('The given conditions cannot be found, these are the',
             'defined conditions: %s')
    txt <- sprintf(paste(txt, collapse = ' '),
                   paste(availableConditions, collapse = ', '))
    stop(txt)
  } else {
    loopCols <- loopCols[match(names(loopCols), conditions)]
  }

  ## get the timepoints for each condition if possible
  timeCol <- .giveMetaoption(x, 'timeCol')
  if (is.null(timeCol)) {
    timepoints.x <- seq_along(loopCols[[1]])
    timepoints.y <- seq_along(loopCols[[2]])
  } else {
    timepoints.x <- colData(x)[loopCols[[1]], timeCol]
    timepoints.y <- colData(x)[loopCols[[2]], timeCol]
  }

  ## if timepoints do not match try to match them
  if (!all(timepoints.x == timepoints.y)) {
    timepoints.x <- timepoints.x[match(timepoints.x, timepoints.y)]
    timepoints.y <- timepoints.y[match(timepoints.x, timepoints.y)]

    mat.x <- mat.x[,which(!is.na(timepoints.x))]
    mat.y <- mat.y[,which(!is.na(timepoints.y))]
    timepoints.x <- timepoints.x[which(!is.na(timepoints.x))]
    timepoints.y <- timepoints.y[which(!is.na(timepoints.y))]
  }

  ## there are no matching timepoints error
  if (length(timepoints.x) == 0 | length(timepoints.y) == 0) {
    txt <- sprintf('The timepoints do not coincide: %s; %s.',
                   paste(timepoints.x, collapse = ', '),
                   paste(timepoints.y, collapse = ', '))
    stop(txt)
  }

  ## make a long format data.frame for plotting
  plotDf <- data.frame(Cond1 = as.vector(mat[, loopCols[[1]]]),
                       Cond2 = as.vector(mat[, loopCols[[2]]]),
                       Time = rep(c(timepoints.x, timepoints.y),
                                  each = nrow(mat)))

  ## remove NAs
  plotDf <- subset(plotDf, !is.na(plotDf$Cond1))
  plotDf <- subset(plotDf, !is.na(plotDf$Cond2))

  ## change column names to conditions
  colnames(plotDf)[seq_len(2)] <- conditions

  if (returnDataFrame) {
    return(plotDf)
  }

  ## plotting ------------------------------------------------------------------

  p <- ggplot(data = plotDf,
              aes_string(x = conditions[1], y = conditions[2])) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'grey70', linetype = 2) +
    facet_wrap(~Time, nrow = 1)  +
    theme_bw()

  p

})

#' @rdname scatterCompareAssays
#' @export
setMethod('scatterCompareAssays', 'PeptideExperiment',
          function(x,
                   conditions,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  callNextMethod()

})

#' @rdname scatterCompareAssays
#' @export
setMethod('scatterCompareAssays', 'ProteomicsExperiment',
          function(x,
                   conditions,
                   assayName,
                   mode = 'protein',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  if (mode == 'protein') {

    scatterCompareAssays(x = x@ProteinExperiment,
                         conditions = conditions,
                         assayName = assayName,
                         returnDataFrame = returnDataFrame,
                         conditionCol = conditionCol,
                         timeCol = timeCol)

  } else if (mode == 'peptide') {

    scatterCompareAssays(x = x@PeptideExperiment,
                         conditions = conditions,
                         assayName = assayName,
                         returnDataFrame = returnDataFrame,
                         conditionCol = conditionCol,
                         timeCol = timeCol)
  }

})
