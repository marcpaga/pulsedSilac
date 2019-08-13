#' @name plotDistributionModel
#' @title Distribution of modelling output
#'
#' @description Plot the distribution of the different model parameters and
#' metrics for each condition.
#'
#' @param modelList A list containing all the model objects, this should be the
#' output of \code{\link{modelTurnover}}.
#' @param value A \code{character} indicating which metric to plot.  Check
#' \code{names(modelList)} for available options. (Default = 'param_values')
#' @param plotType A \code{character} indicating which geometry to plot:
#' 'boxplot' or 'density'. (default = 'density')
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#'
#' @return A \code{ggplot} density or boxplot object, or the \code{data.frame}
#' used to make the plot.
#'
#' @examples
#' data('wormsPE')
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
#' plotDistributionModel(modelList = modelList,
#'                       value = 'param_values',
#'                       plotType = 'density')
#'
#' @export
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
plotDistributionModel <- function(modelList,
                                  value = 'param_values',
                                  plotType = 'density',
                                  returnDataFrame = FALSE) {

  ## argument checker ----------------------------------------------------------
  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  available_metrics <- names(modelList)

  if (!value %in% available_metrics) {
    txt <- c('%s not found in names(modelList), available options are: %s')
    txt <- sprintf(txt, value, paste(available_metrics, collapse = ', '))
    stop(txt)
  }

  if (!plotType %in% c('density', 'boxplot')) {
    txt <- c('%s plotType not valid, please use "density" or "boxplot."')
    txt <- sprintf(txt, plotType)
    stop(txt)
  }

  ## decide type of plot -------------------------------------------------------

  ## there are three categories of plots that can happen here
  ## - parameters and their derivatives (model wide * number of parameters)
  ## - residuals and weights (assay like size)
  ## - error and AIC (model wide)
  n_samples <- length(unlist(attributes(modelList)[['loopCols']]))
  n_conditions <- length(attributes(modelList)[['loopCols']])

  dataToPlot <- modelList[[value]]

  ## parameter type
  if (is.list(dataToPlot)) {
    plotDf <- .plotDistributionModel.Parameter(data = dataToPlot,
                                               ml_attr = attributes(modelList),
                                               returnDataFrame = returnDataFrame,
                                               plotType = plotType,
                                               value = value)
  ## assay like type
  } else if (is.matrix(dataToPlot) & ncol(dataToPlot) == n_samples) {

    plotDf <- .plotDistributionModel.Assay(dataToPlot,
                                           ml_attr = attributes(modelList),
                                           returnDataFrame = returnDataFrame,
                                           plotType = plotType,
                                           value = value)

  ## model wide type
  } else if (is.matrix(dataToPlot) & ncol(dataToPlot) == n_conditions) {

    plotDf <- .plotDistributionModel.Model(dataToPlot,
                                           ml_attr = attributes(modelList),
                                           returnDataFrame = returnDataFrame,
                                           plotType = plotType,
                                           value = value)
  }

  if (returnDataFrame) {
    return(plotDf)
  } else {
    plotDf
  }

}

#' @keywords internal
.plotDistributionModel.Parameter <- function(data,
                                             ml_attr,
                                             returnDataFrame,
                                             plotType,
                                             value) {

  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  value_lab <- value
  for (i in seq_len(length(data))) {
    if (i == 1) {
      dfList <- list()
    }

    cond_vec <- rep(unique(ml_attr[['cond']]),
                    each = nrow(data[[i]]))
    plotDf <- data.frame(value = as.vector(as.matrix(data[[i]])),
                         condition = cond_vec,
                         param = rep(names(data)[i],
                                     times = nrow(data[[i]])))
    dfList[[i]] <- plotDf

  }

  plotDf <- do.call('rbind', dfList)
  plotDf$condition <- as.factor(plotDf$condition)
  plotDf$param <- as.factor(plotDf$param)
  plotDf$condition <- droplevels(plotDf$condition)

  if (returnDataFrame) {
    return(plotDf)
  }

  if (plotType == 'density') {

    p <- ggplot(data = plotDf) +
      geom_density_ridges(aes_string(x = 'value',
                                     y = 'condition',
                                     fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~param, scales = 'free') +
      theme_bw() +
      labs(x = value_lab)

    return(p)

  } else if (plotType == 'boxplot') {

    p <- ggplot(data = plotDf) +
      geom_boxplot(aes_string(x = 'condition',
                              y = 'value',
                              fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~param, scales = 'free') +
      theme_bw() +
      labs(y = value_lab)

    return(p)
  }

}

#' @keywords internal
.plotDistributionModel.Assay <- function(data,
                                         ml_attr,
                                         returnDataFrame,
                                         plotType,
                                         value) {

  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  value_lab <- value
  loopCols <- ml_attr$loopCols
  for (i in seq_len(length(loopCols))) {
    if (i == 1) {
      dfList <- list()
    }

    data_matrix_cond <- data[, loopCols[[i]]]
    time <- ml_attr[['time']][loopCols[[i]]]
    plotDf <- data.frame(value = as.vector(data_matrix_cond),
                         condition = unique(ml_attr[['cond']])[i],
                         time = rep(time,
                                    each = nrow(data_matrix_cond)))

    dfList[[i]] <- plotDf

  }
  plotDf <- do.call('rbind', dfList)
  plotDf$time <- as.factor(plotDf$time)
  plotDf$condition <- as.factor(plotDf$condition)
  plotDf$condition <- droplevels(plotDf$condition)

  if (returnDataFrame) {
    return(plotDf)
  }

  if (plotType == 'density') {

    p <- ggplot(data = plotDf) +
      geom_density_ridges(aes_string(x = 'value',
                                     y = 'time',
                                     fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~condition, scales = 'free') +
      theme_bw() +
      labs(x = value_lab)

    return(p)

  } else if (plotType == 'boxplot') {

    p <- ggplot(data = plotDf) +
      geom_boxplot(aes_string(x = 'time',
                              y = 'value',
                              fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      facet_wrap(~condition, scales = 'free') +
      theme_bw() +
      labs(y = value_lab)

    return(p)
  }

}

#' @keywords internal
.plotDistributionModel.Model <- function(data,
                                         ml_attr,
                                         returnDataFrame,
                                         plotType,
                                         value) {

  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  value_lab <- value
  loopCols <- ml_attr$loopCols
  for (i in seq_len(length(loopCols))) {
    if (i == 1) {
      dfList <- list()
    }

    data_matrix_cond <- data[, i]
    plotDf <- data.frame(value = as.vector(data_matrix_cond),
                         condition = unique(ml_attr[['cond']])[i])

    dfList[[i]] <- plotDf

  }
  plotDf <- do.call('rbind', dfList)
  plotDf$condition <- as.factor(plotDf$condition)
  plotDf$condition <- droplevels(plotDf$condition)

  if (returnDataFrame) {
    return(plotDf)
  }

  plotDf <- subset(plotDf, !is.infinite(value))

  if (plotType == 'density') {

    p <- ggplot(data = plotDf) +
      geom_density_ridges(aes_string(x = 'value',
                                     y = 'condition',
                                     fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      theme_bw() +
      labs(x = value_lab)

    return(p)

  } else if (plotType == 'boxplot') {

    p <- ggplot(data = plotDf) +
      geom_boxplot(aes_string(x = 'condition',
                              y = 'value',
                              fill = 'condition')) +
      scale_fill_manual(values = cbPalette) +
      theme_bw() +
      labs(y = value_lab)

    return(p)
  }


}
