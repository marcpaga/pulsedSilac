#' @rdname scatterCompareModels
#' @name scatterCompareModels
#' @title Scatter plot of two conditions for a model metric.
#'
#' @description Scatter plot of two conditions/replicates for a selected metric
#' of a model. For example to compare turnover rates, model errors...
#  If applicable, timepoints and different paramenters are separated using
#' \code{facet_wrap}.
#'
#' @param modelList A list containing all the model objects, this should be the
#' output of \code{\link{modelTurnover}} with returnModel as TRUE.
#' @param conditions A \code{character} of length 2 indicating which 2
#' conditions should be compared.
#' @param value A \code{character} indicating which metric to plot.  Check
#' \code{names(modelList)} for available options. (Default = 'param_values')
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#'
#' @return A \code{ggplot} object or the \code{data.frame} that would be used
#' instead in the plot.
#'
#' @examples
#' wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
#'
#' modelList <- modelTurnover(x = wormsPE[1:10],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t)',
#'                            start = list(k = 0.02),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' scatterCompareModels(modelList = modelList,
#'                      conditions = c('OW40', 'OW450'),
#'                      value = 'param_values')
#'
#' @import ggplot2
#' @export
scatterCompareModels <- function(modelList,
                                 conditions,
                                 value = 'param_values',
                                 returnDataFrame = FALSE) {

  ## argument checkers ---------------------------------------------------------
  available_metrics <- names(modelList)

  if (!value %in% available_metrics) {
    txt <- c('%s not found in names(modelList), available options are: %s')
    txt <- sprintf(txt, value, paste(available_metrics, collapse = ', '))
    stop(txt)
  }

  ## configuration of metaoption -----------------------------------------------
  if (length(conditions) != 2) {
    stop('conditions must be a character vector of lenght 2')
  }

  loopCols <- attributes(modelList)[['loopCols']]
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


  ## decide type of plot -------------------------------------------------------
  ## there are three categories of plots that can happen here
  ## - parameters and their derivatives (parameter facet_wrap)
  ## - residuals and weights (timepoint facet_wrap)
  ##     - requires timepoint matching
  ## - error and AIC (no facet_wrap)

  n_samples <- length(unlist(attributes(modelList)[['loopCols']]))
  n_conditions <- length(attributes(modelList)[['loopCols']])

  dataToPlot <- modelList[[value]]

  ## parameter type
  if (is.list(dataToPlot)) {

    ## add column names to the data and then subset by column name and order
    ## according to the conditions
    dataToPlot <- lapply(dataToPlot, function(x) {
     colnames(x) <- levels(attributes(modelList)$cond)
     x <- x[, match(colnames(x), conditions)]

     return(x)
    })

    plotDf <- .plotCompareModel.Parameter(dataToPlot,
                                          returnDataFrame)

  ## assay like
  } else if (is.matrix(dataToPlot) & ncol(dataToPlot) == n_samples) {

    ## add column names to the data and then subset by column name and order
    ## according to the conditions
    colnames(dataToPlot) <- attributes(modelList)$cond
    cols <- unlist(lapply(conditions,
                          function(x) which(colnames(dataToPlot) == x)))
    dataToPlot <- dataToPlot[, cols]

    plotDf <- .plotCompareModel.Assay(dataToPlot,
                                      ml_attr = attributes(modelList),
                                      loopCols = loopCols,
                                      returnDataFrame)
  ## model like
  } else if (is.matrix(dataToPlot) & ncol(dataToPlot) == n_conditions) {

    ## add column names to the data and then subset by column name and order
    ## according to the conditions
    colnames(dataToPlot) <- levels(attributes(modelList)$cond)
    dataToPlot <- dataToPlot[, match(colnames(dataToPlot), conditions)]

    plotDf <- .plotCompareModel.Model(dataToPlot,
                                      returnDataFrame)
  }

  if (returnDataFrame) {
    return(plotDf)
  } else {
    plotDf
  }

}


#' @keywords internal
.plotCompareModel.Parameter <- function(data,
                                        returnDataFrame) {

  ## loop over the different parameters
  for (i in seq_len(length(data))) {
    if (i == 1) {
      dfList <- list()
    }

    plotDf <- data.frame(val1 = as.vector(data[[i]][, 1]),
                         val2 = as.vector(data[[i]][, 2]),
                         param = rep(names(data)[i],
                                     times = nrow(data[[i]])))
    dfList[[i]] <- plotDf

  }
  ## join all the parameters into one data.frame
  plotDf <- do.call('rbind', dfList)
  plotDf$param <- as.factor(plotDf$param)
  colnames(plotDf)[1:2] <- colnames(data[[1]])


  if (returnDataFrame) {
    return(plotDf)
  }

  ## remove missing values because they make the coordinate limits look weird
  plotDf <- plotDf[which(apply(plotDf[,1:2], 1,
                               function(x) sum(is.na(x))) < 1),]

  p <- ggplot(data = plotDf) +
    geom_point(aes_string(x = colnames(plotDf)[1],
                          y = colnames(plotDf)[2])) +
    facet_wrap(~param, scales = 'free') +
    geom_abline(slope = 1, intercept = 0, color = 'grey70', linetype = 2) +
    theme_bw()

  return(p)

}

.plotCompareModel.Assay <- function(data, ml_attr, loopCols, returnDataFrame) {

  ## timepoints for matching
  timepoints.x <- ml_attr$time[ml_attr$loopCols[[1]]]
  timepoints.y <- ml_attr$time[ml_attr$loopCols[[2]]]

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
  plotDf <- data.frame(Cond1 = as.vector(data[, loopCols[[1]]]),
                       Cond2 = as.vector(data[, loopCols[[2]]]),
                       Time = rep(c(timepoints.x, timepoints.y),
                                  each = nrow(data)))

  ## remove NAs
  plotDf <- subset(plotDf, !is.na(plotDf$Cond1))
  plotDf <- subset(plotDf, !is.na(plotDf$Cond2))

  ## change column names to conditions
  colnames(plotDf)[1:2] <- names(loopCols_

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

  return(p)

}

#' @keywords internal
.plotCompareModel.Model <- function(data, returnDataFrame) {

  ## no additional processing needed
  plotDf <- as.data.frame(data)

  if (returnDataFrame) {
    return(plotDf)
  }

  ## remove missing values because they make the coordinate limits look weird
  plotDf <- plotDf[which(apply(plotDf[,1:2], 1,
                               function(x) sum(is.na(x))) < 1),]

  p <- ggplot(data = plotDf) +
    geom_point(aes_string(x = colnames(plotDf)[1],
                          y = colnames(plotDf)[2])) +
    geom_abline(slope = 1, intercept = 0, color = 'grey70', linetype = 2) +
    theme_bw()

  return(p)

}
