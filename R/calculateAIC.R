#' Calculates the Akaike Information Criteria (AIC)
#'
#' If only the ProteomicsExperiment argument is given, it will try to find all
#' the computed models and calculate the AIC.
#' If the length of modelName is longer than 1, then all the other paramaters
#' are repeated until the same length of modelName is reached.
#'
#' The following formulas are used to compute the AIC and AICc (small sample
#' size correction):
#' \deqn{AIC = 2k - 2ln(logLik)}
#' \deqn{AICc = AIC + \frac{2k(k + 1)}{n - k - 1}}
#'
#' @param modelList a \code{list} with the model metrics,
#' the output from \code{\link{modelTurnover}}.
#' @param smallSampleSize a \code{logical} indicating if the AIC small
#' sample size correction formula should be used.
#'
#' @return a \code{list} with the model metrics (the given input) plus a matrix
#' named "AIC" with the AIC for each value
#'
#' @seealso \code{\link{compareAIC}},
#'          \code{\link{modelTurnover}}
#' @importFrom stats AIC
#' @importFrom sme AICc
#' @export
calculateAIC <- function(modelList,
                         smallSampleSize = TRUE) {

  if (!is.list(modelList)) {
    stop('"x" must be a list')
  }

  loopCols <- attributes(modelList)[['loopCols']]

  aicMatrix <- matrix(data = NA,
                      nrow = nrow(modelList$stderror),
                      ncol = ncol(modelList$stderror))

  for (i in seq_along(loopCols)) {

    inputList <- modelList$models[[i]]

    ## choose function depending on the sampleSize
    if (smallSampleSize) {
      aic <- lapply(inputList, function(x) {
        if (is.na(x)) {
          NA
        } else {
          AICc(x)
        }
      })
    } else {
      aic <- lapply(inputList, function(x) {
        if (is.na(x)) {
          NA
        } else {
          AIC(x)
        }
      })
    }

    aicMatrix[,i] <- unlist(aic)

  }

  modelList[['AIC']] <- aicMatrix
  return(modelList)

}
