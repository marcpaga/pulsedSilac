#' @rdname modelTurnover
#' @name modelTurnover
#' @title Estimate protein/peptide turnover
#'
#' @description Method to apply turnover models on protein/peptide data
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object.
#' @param assayName \code{character} indicating which assay to use as data
#' input for the model.
#' @param formula \code{formula} to be used. Time must always be named "t" and
#' the data must be named "fraction".
#' @param start named \code{list} with the initical values for the parameters
#' in formula.
#' @param robust \code{logical} indicating if robust modelling from the
#' \code{robustbase} package should be used.
#' @param mode \code{character} indicating which type of data should be used.
#' Can be "protein": one model per protein; "grouped": one model per protein
#' using peptide data; "peptide" one model per peptide.
#' @param verbose \code{logical} indicating if a progress bar should be
#' printed.
#' @param returnModel \code{logical} indicating if the model objects should
#' be returned also in the output.
#' @param conditionCol \code{character} indicating which column of colData(x)
#' describes the conditions.
#' @param timeCol \code{character} indicating which column of colData(x)
#' describes time.
#' @param proteinCol \code{character} indicating which column of rowData(x)
#' describes the assigned protein to a peptide. (Only for peptide data)
#' @param ... further parameters passed into \code{nls} or \code{nlrob}.
#'
#' @return A named \code{list} with either model metrics in matrices or the
#' model objects.
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
#' @importFrom robustbase nlrob
#' @importFrom R.utils insert
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats nls as.formula
#' @import methods
#' @export
setGeneric('modelTurnover', function(x, ...) {
  standardGeneric('modelTurnover')
})

#' @rdname modelTurnover
#' @export
setMethod('modelTurnover',
          'ProteinExperiment',
          function(x,
                   assayName = 'fraction',
                   formula = 'fraction ~ 1-exp(-k*t)',
                   start = list(k = 0.02),
                   robust = FALSE,
                   mode = 'protein',
                   verbose = FALSE,
                   returnModel = FALSE,
                   conditionCol,
                   timeCol,
                   ...){

  ## argument checker ----------------------------------------------------------
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }

  ## data processing and configuration -----------------------------------------

  mat <- assays(x)[[assayName]]

  ## columns of each condition
  loopCols <- .loopWrapper(x, 'conditionCol')
  if (any(vapply(loopCols, length, integer(1)) == 0)) {
    loopCols <- loopCols[which(vapply(loopCols, length, integer(1)) != 0)]
  }

  ## get the condition and time columns to get the vectors from colData
  conditionCol <- .giveMetaoption(x, 'conditionCol')
  timeCol <- .giveMetaoption(x, 'timeCol')
  timeAttr  <- colData(x)[, timeCol]
  condAttr <- colData(x)[, conditionCol]

  ## rownames, colnames and conditionnames for dimension naming of output
  ## matrices
  r_names <- rownames(x)
  c_names <- colnames(x)
  cond_names <- names(loopCols)

  ## if models are returned, then we need a list of lists
  ## else then we need lists and matrices
  if (returnModel) {
    modelList <- list()
  }
  ## initialize all the output matrices ----------------------------------------
  residual_matrix <- matrix(data = NA,
                            nrow = nrow(x),
                            ncol = ncol(x))
  if (!is.null(r_names)) rownames(residual_matrix) <- r_names
  if (!is.null(c_names)) colnames(residual_matrix) <- c_names

  stderror_matrix <- matrix(data = NA,
                            nrow = nrow(x),
                            ncol = length(loopCols))

  if (!is.null(r_names)) rownames(stderror_matrix) <- r_names
  if (!is.null(cond_names)) colnames(stderror_matrix) <- cond_names
  ## there will be a value per condition and parameter, therefore they
  ## are multiplicative
  param_values <- list()
  param_stderror <- list()
  param_tval <- list()
  param_pval <- list()
  for (i in seq_len(length(start))) {
    temp_mat <- matrix(data = NA,
                       nrow = nrow(x),
                       ncol = length(loopCols))

    if (!is.null(r_names)) rownames(temp_mat) <- r_names
    if (!is.null(cond_names)) colnames(temp_mat) <- cond_names

    param_values[[i]] <- temp_mat
    param_stderror[[i]] <- temp_mat
    param_tval[[i]] <- temp_mat
    param_pval[[i]] <- temp_mat
  }
  names(param_values) <- names(start)
  names(param_stderror) <- names(start)
  names(param_tval) <- names(start)
  names(param_pval) <- names(start)


  ## weights are only with robust modelling
  if (robust) {
    weight_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    if (!is.null(r_names)) rownames(weight_matrix) <- r_names
    if (!is.null(c_names)) colnames(weight_matrix) <- c_names
  }

  ## protein turnover modelling ------------------------------------------------
  ## for each condition and for each protein model protein turnover
  for (i in seq_along(loopCols)) {

    if (returnModel) {
      modelList[[i]] <- list()
    }

    if (verbose) {
      cat('Modelling a condition\n')
    }

    for (j in seq_len(nrow(x))) {

      ## progress bar
      if (verbose) {
        if (j == 1){
          pb <- txtProgressBar(min = 1, max = nrow(x), initial = 0, style = 3)
        }
        setTxtProgressBar(pb = pb, value = j)
        if (j == nrow(x)) {
          cat('\n')
        }
      }

      ## modelDf contains the data to do the model on
      if (j == 1) {
        modelDf <- data.frame(t = colData(x)[loopCols[[i]], timeCol],
                              fraction = NA)
      }
      modelDf[, 'fraction'] <- mat[j, loopCols[[i]]]
      modeldata <- .modelTurnover(data = modelDf,
                                  formula = formula,
                                  start = start,
                                  robust = robust,
                                  returnModel = returnModel)

      if (returnModel) {
        if (is.null(modeldata)) {
          modelList[[i]][[j]] <- NA
        } else {
          modelList[[i]][[j]] <- modeldata[['model']]
        }

      }

      if (is.null(modeldata)) {
        next
      }

      ## extract the data from the model and put it in the output matrices
      residual_matrix[j, loopCols[[i]]] <- modeldata[['residuals']]
      stderror_matrix[j, i] <- modeldata[['stderror']]

      for (param in seq_len(length(start))) {
        param_values[[param]][j, i] <- modeldata[['params.vals']][param]
        param_tval[[param]][j, i] <- modeldata[['params.tval']][param]
        param_pval[[param]][j, i] <- modeldata[['params.pval']][param]
        param_stderror[[param]][j, i] <- modeldata[['params.stderror']][param]
      }

      if (robust) {
        weight_matrix[j, loopCols[[i]]] <- modeldata[['weights']]
      }

    } ## end of row loop

    ##residuals and weights as assays with model name
    ## rest in a matrix
  } ## end of loopCols loop


  ## all the output matrices in a list
  outList <- list(residuals = residual_matrix,
                  stderror = stderror_matrix,
                  param_values = param_values,
                  param_pval = param_pval,
                  param_tval = param_tval,
                  param_stderror = param_stderror)

  if (robust) {
    outList[['weights']] <- weight_matrix
  }

  if (returnModel) {
    outList[['models']] <- modelList
  }


  ## add the configuration as attributes that are using in the plotting
  ## functions
  attributes(outList)[['loopCols']] <- loopCols
  attributes(outList)[['time']] <- timeAttr
  attributes(outList)[['cond']] <- condAttr
  attributes(outList)[['assayName']] <- assayName
  attributes(outList)[['mode']] <- mode

  return(outList)
})


#' @rdname modelTurnover
#' @export
setMethod('modelTurnover',
          'PeptideExperiment',
          function(x,
                   assayName = 'fraction',
                   formula = 'fraction ~ 1-exp(-k*t)',
                   start = list(k = 0.02),
                   robust = FALSE,
                   mode = c('grouped', 'peptide'),
                   verbose = FALSE,
                   returnModel = FALSE,
                   conditionCol,
                   timeCol,
                   proteinCol,
                   ...){

  ## argument checker ----------------------------------------------------------
  if (!mode %in% c('grouped', 'peptide')) {
    stop('Mode must be either "grouped" or "peptide".')
  }
  ## a model for each peptide
  if (mode == 'peptide') {
    message('Modelling each peptide individually')
    outList <- callNextMethod()
    return(outList)

  } else if (mode == 'grouped') {
    message('Modelling peptides grouped by protein')
  }

  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  ## metaoptions part
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }
  if (!missing(proteinCol)) {
    metaoptions(x)[['proteinCol']] <- proteinCol
  }

  ## data processing and configuration -----------------------------------------
  mat <- assays(x)[[assayName]]

  ## columns of each condition
  loopCols <- .loopWrapper(x, 'conditionCol')
  if (any(vapply(loopCols, length, integer(1)) == 0)) {
    loopCols <- loopCols[which(vapply(loopCols, length, integer(1)) != 0)]
  }
  ## get the condition and time columns to get the vectors from colData
  conditionCol <- .giveMetaoption(x, 'conditionCol')
  timeCol <- .giveMetaoption(x, 'timeCol')
  proteinCol <- .giveMetaoption(x, 'proteinCol')

  timeAttr  <- colData(x)[, timeCol]
  condAttr <- colData(x)[, conditionCol]
  protAttr <- as.character(rowData(x)[, proteinCol])
  proteinIds <- unique(rowData(x)[, proteinCol])

  ## rownames, colnames and conditionnames for dimension naming of output
  ## matrices
  if (mode == 'peptide') {
    r_names <- rownames(x)
  } else if (mode == 'grouped') {
    ## passed by the ProteomicsExperiment method
    if (exists('r_names_prot')) {
      r_names <- r_names_prot
      r_names_pept <- rownames(x)
    } else {
      r_names <- proteinIds
      r_names_pept <- rownames(x)
    }
  }
  c_names <- colnames(x)
  cond_names <- names(loopCols)

  ## if models are returned, then we need a list of lists
  ## else then we need lists and matrices
  if (returnModel) {
    modelList <- list()
  }

  ## initialize all the output matrices ----------------------------------------
  residual_matrix <- matrix(data = NA,
                            nrow = nrow(x),
                            ncol = ncol(x))
  if (!is.null(r_names_pept)) rownames(residual_matrix) <- r_names_pept
  if (!is.null(c_names)) colnames(residual_matrix) <- c_names

  stderror_matrix <- matrix(data = NA,
                            nrow = length(proteinIds),
                            ncol = length(loopCols))

  if (!is.null(r_names)) rownames(stderror_matrix) <- r_names
  if (!is.null(cond_names)) colnames(stderror_matrix) <- cond_names
  ## there will be a value per condition and parameter, therefore they
  ## are multiplicative
  param_values <- list()
  param_stderror <- list()
  param_tval <- list()
  param_pval <- list()
  for (i in seq_len(length(start))) {

    temp_mat <- matrix(data = NA,
                       nrow = length(proteinIds),
                       ncol = length(loopCols))

    if (!is.null(r_names)) rownames(temp_mat) <- r_names
    if (!is.null(cond_names)) colnames(temp_mat) <- cond_names

    param_values[[i]] <- temp_mat
    param_stderror[[i]] <- temp_mat
    param_tval[[i]] <- temp_mat
    param_pval[[i]] <- temp_mat
  }
  names(param_values) <- names(start)
  names(param_stderror) <- names(start)
  names(param_tval) <- names(start)
  names(param_pval) <- names(start)

  ## weights are only with robust modelling
  if (robust) {
    weight_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    if (!is.null(r_names_pept)) rownames(weight_matrix) <- r_names_pept
    if (!is.null(c_names)) colnames(weight_matrix) <- c_names
  }

  ## protein/peptide turnover modelling ----------------------------------------
  ## for each condition and for each protein model protein turnover
  for (i in seq_along(loopCols)) {

    if (returnModel) {
      modelList[[i]] <- list()
    }

    if (verbose) {
      cat('Modelling a condition\n')
    }

    for (j in seq_along(proteinIds)) {

      ## progress bar
      if (verbose) {
        if (j == 1){
          pb <- txtProgressBar(min = 1, max = length(proteinIds),
                               initial = 0, style = 3)
        }
        setTxtProgressBar(pb = pb, value = j)
        if (j == length(proteinIds)) {
          cat('\n')
        }
      }

      id <- proteinIds[j]
      ## cant use subset because proteinCol is an object
      fracs <- mat[which(rowData(x)[, proteinCol] == id), loopCols[[i]], drop = FALSE]

      ## modelDf contains the data to do the model on
      modelDf <- data.frame(t = rep(colData(x)[loopCols[[i]], timeCol],
                                    each = nrow(fracs)),
                            fraction = as.vector(fracs))

      modeldata <- .modelTurnover(data = modelDf,
                                  formula = formula,
                                  start = start,
                                  robust = robust,
                                  returnModel = returnModel)

      if (returnModel) {
        if (is.null(modeldata)) {
          modelList[[i]][[j]] <- NA
        } else {
          modelList[[i]][[j]] <- modeldata[['model']]
        }

      }

      if (is.null(modeldata)) {
        next
      }

      res <- matrix(modeldata[['residuals']],
                    ncol = length(loopCols[[i]]),
                    nrow = nrow(fracs))

      residual_matrix[which(rowData(x)[,proteinCol] == id), loopCols[[i]]] <- res
      stderror_matrix[j, i] <- modeldata[['stderror']]

      for (param in seq_len(length(start))) {
        param_values[[param]][j, i] <- modeldata[['params.vals']][param]
        param_tval[[param]][j, i] <- modeldata[['params.tval']][param]
        param_pval[[param]][j, i] <- modeldata[['params.pval']][param]
        param_stderror[[param]][j, i] <- modeldata[['params.stderror']][param]
      }

      if (robust) {
        wei <- matrix(modeldata[['weights']],
                      ncol = length(loopCols[[i]]),
                      nrow = nrow(fracs))
        weight_matrix[which(rowData(x)[, proteinCol] == id), loopCols[[i]]] <- wei
      }

    } ## end of row loop

    ##residuals and weights as assays with model name
    ## rest in a matrix
  } ## end of loopCols loop


  ## all the output matrices in a list
  outList <- list(residuals = residual_matrix,
                  stderror = stderror_matrix,
                  param_values = param_values,
                  param_pval = param_pval,
                  param_tval = param_tval,
                  param_stderror = param_stderror)

  if (robust) {
    outList[['weights']] <- weight_matrix
  }

  if (returnModel) {
    outList[['models']] <- modelList
  }

  ## configuration attributes for plotting functions
  attributes(outList)[['loopCols']] <- loopCols
  attributes(outList)[['time']] <- timeAttr
  attributes(outList)[['cond']] <- condAttr
  attributes(outList)[['prot']] <- protAttr
  attributes(outList)[['assayName']] <- assayName
  attributes(outList)[['mode']] <- mode

  return(outList)

})

#' @rdname modelTurnover
#' @export
setMethod('modelTurnover',
          'ProteomicsExperiment',
          function(x,
                   assayName = 'fraction',
                   formula = 'fraction ~ 1-exp(-k*t)',
                   start = list(k = 0.02),
                   robust = FALSE,
                   mode = c('protein', 'grouped', 'peptide'),
                   verbose = FALSE,
                   returnModel = FALSE,
                   conditionCol,
                   timeCol,
                   proteinCol,
                   ...){

  if (!mode %in% c('protein', 'grouped', 'peptide')) {
    stop('Mode must be either "protein", "grouped" or "peptide".')
  }

  if (mode == 'protein') {

    outList <- modelTurnover(x = x@ProteinExperiment,
                             assayName = assayName,
                             formula = formula,
                             start = start,
                             robust = robust,
                             verbose = verbose,
                             returnModel = returnModel,
                             conditionCol = conditionCol,
                             timeCol = timeCol,
                             ...)

  } else if (mode == 'peptide') {

    outList <- modelTurnover(x = x@PeptideExperiment,
                             assayName = assayName,
                             formula = formula,
                             start = start,
                             robust = robust,
                             mode = mode,
                             verbose = verbose,
                             returnModel = returnModel,
                             conditionCol = conditionCol,
                             timeCol = timeCol,
                             proteinCol = proteinCol,
                             ...)

  } else if (mode == 'grouped') {

    outList <- modelTurnover(x = x@PeptideExperiment,
                             assayName = assayName,
                             formula = formula,
                             start = start,
                             robust = robust,
                             mode = mode,
                             verbose = verbose,
                             returnModel = returnModel,
                             conditionCol = conditionCol,
                             timeCol = timeCol,
                             proteinCol = proteinCol,
                             r_names_prot = rownames(x@ProteinExperiment),
                             ...)

  }

  return(outList)
})



#' @keywords internal
.modelTurnover <- function(data, formula, start, robust, returnModel, ...) {

  ## internal function that does the actual modelling, robust or not,
  ## and takes care of NAs

  originalnrow <- nrow(data)

  if (sum(is.na(data[,2])) > 0) {
    isna <- which(!is.na(data[,2]))
    data <- data[isna, ]
  } else {
    isna <- NULL
  }

  if (robust) {
    model  <- try(nlrob(formula = as.formula(formula),
                        data = data,
                        start = start, ...), silent = TRUE)

    if (is(model, 'try-error')) {
      return(NULL)
    }

    residuals2 <- summary(model)[[2]]
    stderror <- summary(model)[[3]]
    weights2 <- summary(model)[[4]]
    params.vals <- summary(model)[[12]][, 1]
    params.stderror <- summary(model)[[12]][, 2]
    params.tval <- summary(model)[[12]][, 3]
    params.pval <- summary(model)[[12]][, 4]

    if (!is.null(isna)) {

      residuals <- rep(NA, originalnrow)
      weights <- rep(NA, originalnrow)
      residuals[isna] <- residuals2
      weights[isna] <- weights2
    } else {
      residuals <- residuals2
      weights <- weights2
    }

    outList <- list(residuals = residuals,
                    stderror = stderror,
                    weights = weights,
                    params.vals = params.vals,
                    params.stderror = params.stderror,
                    params.tval = params.tval,
                    params.pval = params.pval)

    if (returnModel) {
      outList[['model']] <- model
    }

    return(outList)

  } else {
    model  <- try(nls(formula = as.formula(formula),
                      data = data,
                      start = start, ...), silent = TRUE)

    if (is(model, 'try-error')) {
      return(NULL)
    }


    residuals2 <- summary(model)[[2]]
    stderror <- summary(model)[[3]]
    params.vals <- summary(model)[[10]][, 1]
    params.stderror <- summary(model)[[10]][, 2]
    params.tval <- summary(model)[[10]][, 3]
    params.pval <- summary(model)[[10]][, 4]

    if (!is.null(isna)) {
      residuals <- rep(NA, originalnrow)
      residuals[isna] <- residuals2
    }else {
      residuals <- residuals2
    }

    outList <- list(residuals = residuals,
                    stderror = stderror,
                    params.vals = params.vals,
                    params.stderror = params.stderror,
                    params.tval = params.tval,
                    params.pval = params.pval)

    if (returnModel) {
      outList[['model']] <- model
    }
    return(outList)
  }
}
