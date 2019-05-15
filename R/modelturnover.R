#' @export
setGeneric('modelTurnover', function(x, ...) {
  standardGeneric('modelTurnover')
})

#' @rdname modelTurnover
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
#' be the output of the function.
#' @param conditionCol \code{character} indicating which column of colData(x)
#' describes the conditions.
#' @param timeCol \code{character} indicating which column of colData(x)
#' describes time.
#' @param replicateTimeCol \code{character} indicating which column of
#' colData(x) describes the samples time replicates.
#' @param proteinCol \code{character} indicating which column of rowData(x)
#' describes the assigned protein to a peptide. (Only for peptide data)
#' @param ... further parameters passed into \code{nls} or \code{nlrob}.
#' @return A named \code{list} with either model metrics in matrices or the
#' model objects.
#' @importFrom robustbase nlrob
#' @importFrom R.utils insert
#' @importFrom stats nls as.formula
#' @import methods
#' @export
setMethod('modelTurnover',
          'ProteinExperiment',
          function(x,
                   assayName = 'fraction',
                   formula = 'fraction ~ 1-exp(-k*t)',
                   start = list(k = 0.02),
                   robust = FALSE,
                   verbose = FALSE,
                   returnModel = FALSE,
                   conditionCol,
                   timeCol,
                   replicateTimeCol,
                   ...){

  ## argument check
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  mat <- assays(x)[[assayName]]

  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }

  ## which columns belong to which experiment
  loopList <- tryCatch(
    {
      loop <- experimentLoopWrapper(x, 'cond.timerep')
      conditionCol <- giveMetaoption(x, 'conditionCol')
      timeRepCol <- giveMetaoption(x, 'replicateTimeCol')
      cond <- paste(as.character(colData(x)[,conditionCol]),
                    as.character(colData(x)[,timeRepCol]), sep = '.')
      list(loop, cond)
    },
    error = function(c){
      tryCatch(
        {
          loop <- experimentLoopWrapper(x, 'cond')
          conditionCol <- giveMetaoption(x, 'conditionCol')
          cond <- as.character(colData(x)[,conditionCol])
          list(loop, cond)
        },
        error = function(c){
          list(seq_len(ncol(x)), NA)
        }
      )
    }
  )

  ## where to get the timepoint data from
  timeCol <- giveMetaoption(x, 'timeCol')

  ## attributes for plotting
  loopCols <- loopList[[1]]
  ## conditionNames
  condAttr <- loopList[[2]]
  ## time
  timeAttr  <- colData(x)[,timeCol]

  ## if models are returned, then we need a list of lists
  ## else then we need lists and matrices
  if (returnModel) {
    modelList <- list()
  } else {

    residual_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    stderror_matrix <- matrix(data = NA,
                              nrow = nrow(x),
                              ncol = length(loopCols))
    ## there will be a value per condition and parameter, therefore they
    ## are multiplicative
    param_values <- matrix(data = NA,
                           nrow = nrow(x),
                           ncol = length(loopCols)*length(start))
    param_stderror <- matrix(data = NA,
                             nrow = nrow(x),
                             ncol = length(loopCols)*length(start))
    param_tval <- matrix(data = NA,
                         nrow = nrow(x),
                         ncol = length(loopCols)*length(start))
    param_pval <- matrix(data = NA,
                         nrow = nrow(x),
                         ncol = length(loopCols)*length(start))
    ## weights are only with robust modelling
    if (robust) {
      weight_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    }

  }


  ## for each condition and for each protein model protein turnover
  for (i in seq_along(loopCols)) {

    if (returnModel) {
      modelList[[i]] <- list()
    }

    if (verbose) {
      cat('Modelling a condition\n')
    }

    for (j in seq_len(nrow(x))) {

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
          modelList[[i]][[j]] <- modeldata
        }

      } else {

        if (is.null(modeldata)) {
          next
        }

        chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

        residual_matrix[j, loopCols[[i]]] <- modeldata[['residuals']]
        stderror_matrix[j, i] <- modeldata[['stderror']]

        chunks <- chunk(seq_len(ncol(param_values)), length(loopCols))
        param_values[j, chunks[[i]]] <- modeldata[['params.vals']]
        param_tval[j, chunks[[i]]] <- modeldata[['params.tval']]
        param_pval[j, chunks[[i]]] <- modeldata[['params.pval']]
        param_stderror[j, chunks[[i]]] <- modeldata[['params.stderror']]

        if (robust) {
          weight_matrix[j, loopCols[[i]]] <- modeldata[['weights']]
        }

      }

    } ## end of row loop

    ##residuals and weights as assays with model name
    ## rest in a matrix
  } ## end of loopCols loop

  if (returnModel) {
    attributes(modelList)[['loopCols']] <- loopCols
    attributes(modelList)[['time']] <- timeAttr
    attributes(modelList)[['cond']] <- condAttr
    return(modelList)
  }

  outList <- list(residuals = residual_matrix,
                  stderror = stderror_matrix,
                  param_values = param_values,
                  param_pval = param_pval,
                  param_tval = param_tval,
                  param_stderror = param_stderror)

  if (robust) {
    outList[['weights']] <- weight_matrix
  }


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
                   replicateTimeCol,
                   proteinCol,
                   ...){

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

  mat <- assays(x)[[assayName]]

  ## metaoptions part
  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(timeCol)) {
    metaoptions(x)[['timeCol']] <- timeCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }
  if (!missing(proteinCol)) {
    metaoptions(x)[['proteinCol']] <- proteinCol
  }

  ## which columns belong to which experiment
  loopList <- tryCatch(
    {
      loop <- experimentLoopWrapper(x, 'cond.timerep')
      conditionCol <- giveMetaoption(x, 'conditionCol')
      timeRepCol <- giveMetaoption(x, 'replicateTimeCol')
      cond <- paste(as.character(colData(x)[,conditionCol]),
                    as.character(colData(x)[,timeRepCol]), sep = '.')
      list(loop, cond)
    },
    error = function(c){
      tryCatch(
        {
          loop <- experimentLoopWrapper(x, 'cond')
          conditionCol <- giveMetaoption(x, 'conditionCol')
          cond <- as.character(colData(x)[,conditionCol])
          list(loop, cond)
        },
        error = function(c){
          list(seq_len(ncol(x)), NA)
        }
      )
    }
  )

  ## where to get the timepoint data from
  timeCol <- giveMetaoption(x, 'timeCol')
  ## group peptides from the same protein and return a model
  ## this requires the proteinCol metaoption
  proteinCol <- giveMetaoption(x, 'proteinCol')
  proteinIds <- unique(rowData(x)[, proteinCol])

  ## attributes for plotting
  loopCols <- loopList[[1]]
  ## conditionNames
  condAttr <- loopList[[2]]
  ## time
  timeAttr  <- colData(x)[, timeCol]
  ## protein
  protAttr <- as.character(rowData(x)[, proteinCol])

  ## if models are returned, then we need a list of lists
  ## else then we need lists and matrices
  if (returnModel) {
    modelList <- list()
  } else {

    residual_matrix <- matrix(data = NA,
                              nrow = nrow(x),
                              ncol = ncol(x))
    stderror_matrix <- matrix(data = NA,
                              nrow = length(proteinIds),
                              ncol = length(loopCols))
    ## there will be a value per condition and parameter, therefore they
    ## are multiplicative
    param_values <- matrix(data = NA,
                           nrow = length(proteinIds),
                           ncol = length(loopCols)*length(start))
    param_stderror <- matrix(data = NA,
                             nrow = length(proteinIds),
                             ncol = length(loopCols)*length(start))
    param_tval <- matrix(data = NA,
                         nrow = length(proteinIds),
                         ncol = length(loopCols)*length(start))
    param_pval <- matrix(data = NA,
                         nrow = length(proteinIds),
                         ncol = length(loopCols)*length(start))
    ## weights are only with robust modelling
    if (robust) {
      weight_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    }

  }


  ## for each condition and for each protein model protein turnover
  for (i in seq_along(loopCols)) {

    if (returnModel) {
      modelList[[i]] <- list()
    }

    if (verbose) {
      cat('Modelling a condition\n')
    }

    for (j in seq_along(proteinIds)) {

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
      fracs <- mat[which(rowData(x)[,proteinCol] == id), loopCols[[i]], drop = FALSE]

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
          modelList[[i]][[j]] <- modeldata
        }

      } else {

        if (is.null(modeldata)) {
          next
        }

        chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

        res <- matrix(modeldata[['residuals']],
                      ncol = length(loopCols[[i]]),
                      nrow = nrow(fracs))

        residual_matrix[which(rowData(x)[,proteinCol] == id), loopCols[[i]]] <- res
        stderror_matrix[j, i] <- modeldata[['stderror']]

        chunks <- chunk(seq_len(ncol(param_values)), length(loopCols))
        param_values[j, chunks[[i]]] <- modeldata[['params.vals']]
        param_tval[j, chunks[[i]]] <- modeldata[['params.tval']]
        param_pval[j, chunks[[i]]] <- modeldata[['params.pval']]
        param_stderror[j, chunks[[i]]] <- modeldata[['params.stderror']]

        if (robust) {
          wei <- matrix(modeldata[['weights']],
                        ncol = length(loopCols[[i]]),
                        nrow = nrow(fracs))
          weight_matrix[which(rowData(x)[,proteinCol] == id), loopCols[[i]]] <- wei
        }

      }

    } ## end of row loop

    ##residuals and weights as assays with model name
    ## rest in a matrix
  } ## end of loopCols loop

  if (returnModel) {
    attributes(modelList)[['loopCols']] <- loopCols
    attributes(modelList)[['time']] <- timeAttr
    attributes(modelList)[['cond']] <- condAttr
    attributes(modelList)[['prot']] <- protAttr
    return(modelList)
  }

  outList <- list(residuals = residual_matrix,
                  stderror = stderror_matrix,
                  param_values = param_values,
                  param_pval = param_pval,
                  param_tval = param_tval,
                  param_stderror = param_stderror)

  if (robust) {
    outList[['weights']] <- weight_matrix
  }

  attributes(outList)[['loopCols']] <- loopCols
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
                   mode = c('protein','grouped', 'peptide'),
                   verbose = FALSE,
                   returnModel = FALSE,
                   conditionCol,
                   timeCol,
                   replicateTimeCol,
                   proteinCol,
                   ...){

  if (!mode %in% c('protein', 'grouped', 'peptide')) {
    stop('Mode must be either "protein" ,"grouped" or "peptide".')
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
                             replicateTimeCol = replicateTimeCol,
                             ...)

  } else if (mode %in% c('grouped', 'peptide')) {

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
                             replicateTimeCol = replicateTimeCol,
                             proteinCol = proteinCol,
                             ...)

  }

  return(outList)
})



#' @keywords internal
.modelTurnover <- function(data, formula, start, robust, returnModel, ...) {

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
    if (returnModel) {
      return(model)
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

    return(outList)

  } else {
    model  <- try(nls(formula = as.formula(formula),
                      data = data,
                      start = start, ...), silent = TRUE)

    if (is(model, 'try-error')) {
      return(NULL)
    }

    if (returnModel) {
      return(model)
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

    return(outList)
  }
}
