#' @rdname compareAssay
#' @name compareAssay
#' @title Scatter plot of two conditions for each timepoint of an assay.
#'
#' @description Scatter plot of two conditions/replicates for a selected assay.
#' Timepoints are separated using \code{facet_wrap}.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param y A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
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
#' compareAssay(x = wormsPE[, 1:7],
#'              y = wormsPE[, 8:14],
#'              assayName = 'ratio',
#'              mode = 'protein')
#'
#' @import ggplot2
#' @export
setGeneric('compareAssay', function(x, ...){
  standardGeneric('compareAssay')
})

#' @rdname compareAssay
#' @export
setMethod('compareAssay', 'ProteinExperiment',
          function(x,
                   y,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  if (!assayName %in% names(assays(y))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## count how many proteins per sample
  mat.x <- assays(x)[[assayName]]
  mat.y <- assays(y)[[assayName]]

  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
    metaoptions(y)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
    metaoptions(y)[['timeCol']] <- timeCol
  }

  ## define the timepoints
  if (!is.na(metaoptions(x)[['timeCol']])) {
    timepoints.x <- colData(x)[, metaoptions(x)[['timeCol']]]
  } else {
    timepoints.x <- seq_len(nrow(colData(x)))
  }
  if (!is.na(metaoptions(y)[['timeCol']])) {
    timepoints.y <- colData(y)[, metaoptions(y)[['timeCol']]]
  } else {
    timepoints.y <- seq_len(nrow(colData(y)))
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

    timepoints.x <- colData(x)[,metaoptions(x)[['timeCol']]]
    timepoints.y <- colData(y)[,metaoptions(y)[['timeCol']]]

    txt <- sprintf('The timepoints do not coincide: %s; %s.',
                   paste(timepoints.x, collapse = ', '),
                   paste(timepoints.y, collapse = ', '))
    stop(txt)
  }

  ## make a long format data.frame for plotting
  for (t in seq_along(timepoints.x)) {
    if (t == 1) {
      tempList <- list()
    }
    tempDf <- data.frame(Cond1 = as.vector(mat.x[,t]),
                         Cond2 = as.vector(mat.y[,t]),
                         Time = timepoints.x[t])
    tempList[[t]] <- tempDf
  }
  plotDf <- do.call('rbind', tempList)

  ## remove NAs
  plotDf <- subset(plotDf, !is.na(plotDf$Cond1))
  plotDf <- subset(plotDf, !is.na(plotDf$Cond2))

  ## if possible change the column names to the conditions
  if (!is.na(metaoptions(x)[['conditionCol']])) {
    conditionName1 <- unique(colData(x)[,metaoptions(x)[['conditionCol']]])
    conditionName1 <- as.character(conditionName1)
    colnames(plotDf)[1] <- conditionName1
  } else {
    conditionName1 <- 'Cond1'
  }

  if (!is.na(metaoptions(y)[['conditionCol']])) {
    conditionName2 <- unique(colData(y)[,metaoptions(y)[['conditionCol']]])
    conditionName2 <- as.character(conditionName2)
    colnames(plotDf)[2] <- conditionName2
  } else {
    conditionName2 <- 'Cond2'
  }

  if (returnDataFrame) {
    return(plotDf)
  }

  p <- ggplot(data = plotDf,
              aes_string(x = conditionName1, y = conditionName2)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, color = 'grey70', linetype = 2) +
    facet_wrap(~Time, nrow = 1)  +
    theme_bw()

  p

})

#' @rdname compareAssay
#' @export
setMethod('compareAssay', 'PeptideExperiment',
          function(x,
                   y,
                   assayName,
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  callNextMethod()

})

#' @rdname compareAssay
#' @export
setMethod('compareAssay', 'ProteomicsExperiment',
          function(x,
                   y,
                   assayName,
                   mode = 'protein',
                   returnDataFrame = FALSE,
                   conditionCol,
                   timeCol) {

  if (mode == 'protein') {

    compareAssay(x = x@ProteinExperiment,
                 y = y@ProteinExperiment,
                 assayName = assayName,
                 returnDataFrame = returnDataFrame,
                 conditionCol = conditionCol,
                 timeCol = timeCol)

  } else if (mode == 'peptide') {

    compareAssay(x = x@PeptideExperiment,
                 y = y@PeptideExperiment,
                 assayName = assayName,
                 returnDataFrame = returnDataFrame,
                 conditionCol = conditionCol,
                 timeCol = timeCol)
  }

})
