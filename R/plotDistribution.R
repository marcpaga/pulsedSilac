#' @export
setGeneric('plotDistribution', function(x, ...){
  standardGeneric('plotDistribution')
})


#' @title Distribution of modelling output
#'
#' @description Plot the distribution of the different model parameters and
#' metrics for each condition.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object.
#' @param modelList A list containing all the model objects, this should be the
#' output of \link{\code{modelTurnover}} with returnModel as TRUE.
#' @param value A \code{character} indicating which metric to plot: 'parameter',
#' 'error', 'residuals' or 'weights'. (default = 'parameter')
#' @param plotType A \code{character} indicating which geometry to plot:
#' 'boxplot' or 'density'. (default = 'density')
#' @param return A \code{character} indicating what to return a 'plot'
#' or 'data.frame'. (default = 'plot')
#'
#' @return A scatter plot with a fitted line or a \code{data.frame}.
#' @export
#' @import ggplot2
#' @importFrom ggridges geom_density_ridges
setMethod('plotDistribution',
          'ProteinExperiment',
          function(x,
                   modelList,
                   value = 'parameter',
                   plotType = 'density',
                   return = 'plot') {

  ## cb palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
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

    if (return == 'data.frame') {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes(x = value, y = condition, fill = condition)) +
        theme_classic() +
        facet_wrap(~param) +
        scale_fill_manual(values = cbPalette)

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes(x = condition, y = value, fill = condition)) +
        theme_classic()+
        facet_wrap(~param) +
        scale_fill_manual(values = cbPalette)

    }


  }

  ## to plot model errors
  if (value == 'error') {

    cond_vec <- rep(unique(attributes(modelList)[['cond']]),
                    each = nrow(modelList[['stderror']]))
    plotDf <- data.frame(value = as.vector(modelList[['stderror']]),
                         condition = cond_vec)
    plotDf$condition <- as.factor(plotDf$condition)

    if (return == 'data.frame') {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes(x = value, y = condition, fill = condition)) +
        theme_classic() +
        scale_fill_manual(values = cbPalette)

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes(x = condition, y = value, fill = condition)) +
        theme_classic() +
        scale_fill_manual(values = cbPalette)

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

    if (return == 'data.frame') {
      return(plotDf)
    }

    if (plotType == 'density') {

      p <- ggplot(data = plotDf) +
        geom_density_ridges(aes(x = value,
                                y = time,
                                fill = condition)) +
        theme_classic() +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette)

    }

    if (plotType == 'boxplot') {

      p <- ggplot(data = plotDf) +
        geom_boxplot(aes(x = time,
                         y = value,
                         fill = condition)) +
        theme_classic() +
        facet_wrap(~condition) +
        scale_fill_manual(values = cbPalette)

    }

  }

  p

})

#' @export
setMethod('plotDistribution',
          'PeptideExperiment',
          function(x,
                   modelList,
                   value = 'parameter',
                   plotType = 'density',
                   return = 'plot') {


})


#' @export
setMethod('plotDistribution',
          'ProteomicsExperiment',
          function(x,
                   modelList,
                   value = 'parameter',
                   plotType = 'density',
                   return = 'plot') {

  if (attributes(modelList)[['mode']] == 'protein') {

    plotDistribution(x = x@ProteinExperiment,
                     modelList = modelList,
                     value = value,
                     plotType = plotType,
                     return = return)

  } else if (attributes(modelList)[['mode']] == 'peptide') {

    plotDistribution(x = x@PeptideExperiment,
                     modelList = modelList,
                     value = value,
                     plotType = plotType,
                     return = return)

  } else if (attributes(modelList)[['mode']] == 'grouped') {

    plotDistribution(x = x@PeptideExperiment,
                     modelList = modelList,
                     value = value,
                     plotType = plotType,
                     return = return)

  }

})


