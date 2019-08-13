#' @title Merge several models lists.
#'
#' @description Merges several models lists into one. The models lists must come
#' from \code{\link{modelTurnover}} with the same arguments with the exception
#' of the input data. This function should be used in cases where a model is
#' fit for different conditions of an experiment with NO overlap in the samples
#' and NO missing samples. Otherwise the plotting functions might give
#' incorrect outputs.
#'
#' @details When merging the attributes are also merged. Some of these need
#' to be recalculated since they contain information about the input data
#' positioning (column number). These attributes are used in the plotting
#' functions.
#'
#' Take this into consideration so that the order of the models lists follows
#' the order of the columns in the original data. This also means no skipped
#' conditions and no skipped samples. If that is the case, build an
#' intermediari object that contains only the samples to be used (see examples).
#'
#' @param ... Lists with model data, output from \code{\link{modelTurnover}}.
#'
#' @return A list of models.
#'
#' @examples
#' data('wormsPE')
#' wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
#'
#' modelList1 <- modelTurnover(x = wormsPE[1:10, 1:7],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t)',
#'                            start = list(k = 0.02),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' modelList2 <- modelTurnover(x = wormsPE[1:10, 8:14],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t)',
#'                            start = list(k = 0.02),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' mergedModelList <- mergeModelsLists(modelList1, modelList2)
#'
#' @seealso \code{\link{modelTurnover}}
#' @export
mergeModelsLists <- function(...) {

  inputList <- list(...)

  ## input checker -------------------------------------------------------------
  ## need minimum two lists to merge them
  if (length(inputList) < 2) {
    txt <- sprintf(paste0('To merge models lists at least 2 lists are',
                          ' required, %i given.'),
                   length(inputList))
    stop(txt)
  }

  ## check that all the lists have the same names
  namesList <- lapply(inputList, names)
  if (!all(vapply(namesList, FUN.VALUE = logical(1),
                  FUN = identical, namesList[[1]]))) {
    txt <- sprintf('Not all the lists have the same names')
    stop(txt)
  }

  ## check that all the lists have the same matrix nrow
  nrowList <- lapply(inputList, function(x) nrow(x[['residuals']]))
  if (!all(vapply(nrowList, FUN.VALUE = logical(1),
                  FUN = identical, nrowList[[1]]))) {
    txt <- sprintf('Not all the lists have the number of models')
    stop(txt)
  }

  ## check that all the lists have the same parameters
  paramaterList <- lapply(inputList, function(x) names(x[['param_values']]))
  if (!all(vapply(paramaterList, FUN.VALUE = logical(1),
                  FUN = identical, paramaterList[[1]]))) {
    txt <- sprintf('Not all the lists have the same parameters')
    stop(txt)
  }

  ## check that all the lists the same mode
  modeList <- lapply(inputList, function(x) attributes(x)[['mode']])
  if (!all(vapply(modeList, FUN.VALUE = logical(1),
                  FUN = identical, modeList[[1]]))) {
    txt <- sprintf(paste0('Not all the lists have the same mode ("protein", ',
                          '"peptide", or "grouped")'))
    stop(txt)
  }

  ## data processing -----------------------------------------------------------
  ## the first model list will be used as template
  outList <- inputList[[1]]

  for (n in names(outList)) {

    ## if it is a matrix this is enough
    if (is.matrix(outList[[n]])) {
      outList[[n]] <- do.call('cbind', lapply(inputList, '[[', n))

    ## returnModel is TRUE
    } else if (n == 'models') {

      outList[[n]] <- lapply(inputList,
                             FUN = function(x) unlist(x[[n]],
                                                      recursive = FALSE))

    ## otherwise it is a paramater list and it has to be done per paramater
    } else {
      for (param in names(outList[[n]])) {
        m <- do.call('cbind', lapply(inputList, function(x) x[[n]][[param]]))
        outList[[n]][[param]] <- m
      }
    }
  }

  ## finally mix the attributes, from the possible attributes:
  ## - names: no need to change
  ## - loopCols: add to list plus change value
  ## - time: just concatenate
  ## - cond : just concatenate
  ## - assayName: no need to change
  ## - mode: no need to change
  ## - prot: no need to change

  ## concatenate these two attributes
  timeAttr <- lapply(inputList, function(x) attributes(x)[['time']])
  attributes(outList)[['time']] <- do.call('c', timeAttr)

  condAttr <- lapply(inputList, function(x) attributes(x)[['cond']])
  attributes(outList)[['cond']] <- unlist(condAttr)

  ## for loop for the loopCols attribute
  for (i in seq_along(inputList)) {
    if (i == 1) {
      loopCols <- attributes(inputList[[i]])[['loopCols']]
      cols <- length(loopCols[[i]])
      next
    }


    loopCols[[i]] <- unname(unlist(attributes(inputList[[i]])[['loopCols']]))
    loopCols[[i]] <- loopCols[[i]] + cols
    names(loopCols)[i] <- names(attributes(inputList[[i]])[['loopCols']])
  }
  attributes(outList)[['loopCols']] <- loopCols

  return(outList)
}
