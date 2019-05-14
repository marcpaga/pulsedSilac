setGeneric('modelTurnover', function(x, ...) {
  standardGeneric('modelTurnover')
})

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
  loopCols <- tryCatch(
    {
      experimentLoopWrapper(x, 'cond.timerep')
    },
    error = function(c){
      tryCatch(
        {
          experimentLoopWrapper(x, 'cond')
        },
        error = function(c){
          list(seq_len(ncol(x)))
        }
      )
    }
  )
  ## where to get the timepoint data from
  timeCol <- giveMetaoption(x, 'timeCol')

  ## if models are returned, then we need a list of lists
  ## else then we need lists and matrices
  if (returnModel) {
    modelList <- list()
  } else {

    residual_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
    stderror_matrix <- matrix(data = NA, nrow = nrow(x), ncol = ncol(x))
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

    if (robust) {
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
          modeList[[i]][[j]] <- NA
        } else {
          modeList[[i]][[j]] <- modeldata
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





#' Calculates the exponential decay constant
#'
#' @param data dataframe containing time and fraction information
#' @param robust boolean, (default = TRUE)
#' @param formula formula that will be used in the model, note that time must
#' be called 't' and fraction must be called 'fraction',
#' (default = 'fraction ~ 1-exp(-k*t)')
#' @param start initial values for the paramaters to be tuned in the model,
#' (default = list(k = 0.02))
#' @param returnModel logical, if TRUE then the model object is returned
#' @param ... further parameters passed into nls or nlrob
#' @import robustbase
#' @import methods
#' @importFrom R.utils insert
#' @importFrom stats nls as.formula
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
