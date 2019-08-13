#' Calculates the Akaike Information Criteria (AIC)
#'
#' Calculates the AIC for each of the computed models. Requires that
#' \code{modelTurnover} is run with \code{reuturnModel = TRUE}.
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
#' modelList <- calculateAIC(modelList, smallSampleSize = TRUE)
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

  if (!'models' %in% names(modelList)) {
    stop('There are no models in this modelList, did you run modelTurnover ',
         'with returnModel = TRUE?')
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
        if (is(x, 'nls')) {
          AICc(x)
        } else {
          NA
        }
      })
    } else {
      aic <- lapply(inputList, function(x) {
        if (is(x, 'nls')) {
          AIC(x)
        } else {
          NA
        }
      })
    }

    aicMatrix[,i] <- unlist(aic)

  }

  modelList[['AIC']] <- aicMatrix
  return(modelList)

}
