#' @title Distribution of modelling output
#'
#' @description Plot the distribution of the different model parameters and
#' metrics for each condition.
#'
#' @param modelList A list containing all the model objects, this should be the
#' output of \code{\link{modelTurnover}} with returnModel as TRUE.
#' @param value A \code{character} indicating which metric to plot: 'parameter',
#' 'error', 'residuals', 'weights', 'aicvalue' or 'aicprobabilities'.
#' (default = 'parameter')
#' @param plotType A \code{character} indicating which geometry to plot:
#' 'boxplot' or 'density'. (default = 'density')
#' @param returnDataFrame A \code{logical} indicating if the \code{data.frame}
#' used for the plot should be returned instead.
#'
#' @return A \code{ggplot} density or boxplot object, or the \code{data.frame}
#' used to make the plot.
#'
#' @examples
#'
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
#' plotDistribution(modelList = modelList,
#'                  value = 'error',
#'                  plotType = 'density')
#'
#' @export
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
plotDistribution <- function(modelList,
                             value = 'parameter',
                             plotType = 'density',
                             returnDataFrame = FALSE) {

  ## cb palette
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

  ## to plot the different parameters
  if (value == 'parameter') {

    for (i in seq_len(length(modelList[['param_values']]))) {
      if (i == 1) {
        dfList <- list()
      }
      param_list <- modelList[['param_values']]

      cond_vec <- rep(unique(attributes(modelList)[['cond']]),
                      each = nrow(param_list[[i]]))
      plotDf <- data.frame(value = as.vector(param_list[[i]]),
                           condition = cond_vec,
                           param = rep(names(param_list)[i],
                                       times = nrow(param_list[[i]])))
      dfList[[i]] <- plotDf

    }
    plotDf <- do.call('rbind', dfList)
    plotDf$condition <- as.factor(plotDf$condition)
    plotDf$param <- as.factor(plotDf$param)

    if (returnDataFrame) {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes_string(x = 'value',
                                       y = 'condition',
                                       fill = 'condition')) +
        facet_wrap(~param) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes_string(x = 'condition',
                                y = 'value',
                                fill = 'condition')) +
        facet_wrap(~param) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }


  }

  ## to plot model errors
  if (value %in% c('error', 'aicvalue')) {

    if (value == 'error') {
      cond_vec <- rep(unique(attributes(modelList)[['cond']]),
                      each = nrow(modelList[['stderror']]))
      plotDf <- data.frame(value = as.vector(modelList[['stderror']]),
                           condition = cond_vec)
      plotDf$condition <- as.factor(plotDf$condition)
    } else if (value == 'aicvalue') {
      cond_vec <- rep(unique(attributes(modelList)[['cond']]),
                      each = nrow(modelList[['AIC']]))
      plotDf <- data.frame(value = as.vector(modelList[['AIC']]),
                           condition = cond_vec)
      plotDf$condition <- as.factor(plotDf$condition)
    }


    if (returnDataFrame) {
      return(plotDf)
    }

    if (sum(plotDf$value == Inf, na.rm = TRUE) > 0) {
      plotDf$value[plotDf$value == Inf] <- NA
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes_string(x = 'value',
                                       y = 'condition',
                                       fill = 'condition')) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes_string(x = 'condition',
                                y = 'value',
                                fill = 'condition')) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }

  }

  ## to plot the residuals
  if (value %in% c('residuals', 'weights')) {

    if (value == 'residuals') {
      data_matrix <- modelList[['residuals']]
    } else if (value == 'weights') {
      data_matrix <- modelList[['weights']]
    }

    loopCols <- attributes(modelList)$loopCols
    for (i in seq_len(length(loopCols))) {
      if (i == 1) {
        dfList <- list()
      }

      data_matrix_cond <- data_matrix[, loopCols[[i]]]
      time <- attributes(modelList)[['time']][loopCols[[i]]]
      plotDf <- data.frame(value = as.vector(data_matrix_cond),
                           condition = unique(attributes(modelList)[['cond']])[i],
                           time = rep(time,
                                      each = nrow(data_matrix_cond)))

      dfList[[i]] <- plotDf

    }
    plotDf <- do.call('rbind', dfList)
    plotDf$time <- as.factor(plotDf$time)
    plotDf$condition <- as.factor(plotDf$condition)

    if (returnDataFrame) {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes_string(x = 'value',
                                       y = 'time',
                                       fill = 'condition')) +
        theme_classic() +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette)

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes_string(x = 'time',
                                y = 'value',
                                fill = 'condition')) +
        theme_classic() +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette)

    }
  }


  if (value == 'aicprobabilities') {

    n_models <- ncol(modelList[[1]])
    loopCols <- attributes(modelList)[['loopCols']]

    model_names <- paste0('model', seq_len(n_models))

    for (i in seq_len(length(loopCols))) {
      if (i == 1) {
        dfList <- list()
      }

      aic_matrix <- modelList[[i]]
      plotDf <- data.frame(value = as.vector(aic_matrix),
                           model = rep(model_names, each = nrow(aic_matrix)))
      plotDf$condition <- unique(attributes(modelList)[['cond']])[i]

      dfList[[i]] <- plotDf

    }

    plotDf <- do.call('rbind', dfList)
    plotDf$model <- as.factor(plotDf$model)
    plotDf$condition <- as.factor(plotDf$condition)

    if (returnDataFrame) {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes_string(x = 'value',
                                       y = 'model',
                                       fill = 'condition')) +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes_string(x = 'time',
                                y = 'model',
                                fill = 'condition')) +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette) +
        theme_bw()

    }
  }

  p

}

