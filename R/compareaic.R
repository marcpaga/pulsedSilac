#' Calculates the probability of each model from a set of models.
#'
#' For a given set of AIC from models, the probability of each model relative
#' to the rest of the models of the set is calculated using the following
#' formula:
#' \deqn{\prod AIC_{i} = \frac{exp(\frac{AIC_{min}-AIC_{i}}{2})}{\sum_{j}exp(\frac{AIC_{min}-AIC_{j}}{2})}}
#'
#' @param ... a \code{list} with the model metrics,
#' the output from \code{\link{modelTurnover}} and \code{\link{calculateAIC}}.
#' @return a \code{list} with a matrix for each experiment condition. The matrix
#' contains the probabilities of each model (columns) for each protein/peptide
#' (rows).
#'
#' @examples
#'
#' wormsPE <- calculateIsotopeFraction(wormsPE, ratioAssay = 'ratio')
#'
#' modelList1 <- modelTurnover(x = wormsPE[1:10],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t)',
#'                            start = list(k = 0.02),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' modelList1 <- calculateAIC(modelList1, smallSampleSize = TRUE)
#'
#' modelList2 <- modelTurnover(x = wormsPE[1:10],
#'                            assayName = 'fraction',
#'                            formula = 'fraction ~ 1 - exp(-k*t) + b',
#'                            start = list(k = 0.02, b = 0),
#'                            mode = 'protein',
#'                            robust = FALSE,
#'                            returnModel = TRUE)
#'
#' modelList2 <- calculateAIC(modelList2, smallSampleSize = TRUE)
#'
#' modelProbabilities <- compareAIC(modelList1, modelList2)
#'
#' @seealso \code{\link{calculateAIC}},
#'          \code{\link{modelTurnover}}.
#' @export
compareAIC <- function(...) {

  inputLists <- list(...)
  n_models <- length(inputLists)

  aicsList <- lapply(inputLists, '[[', 'AIC')

  n_conditions <- length(attributes(inputLists[[1]])$loopCols)
  outputList <- list()

  for (i in seq_len(n_conditions)) {
    for (j in seq_len(n_models)) {
      if (j == 1) {
        tempList <- list()
      }
      tempList[[j]] <- aicsList[[j]][,i]
    }
    aicMatrix <- do.call('cbind', tempList)
    probs <- t(apply(aicMatrix, 1, .compareAIC))
    outputList[[i]] <- probs
  }

  attrList <- attributes(inputLists[[1]])
  attrList <- attrList[-which(names(attrList) == 'names')]
  attributes(outputList) <- attrList
  return(outputList)
}

#' Calculates the probability of a set of AIC values
#'
#' @param aics Numeric vector containing the AIC values
#' @return A numeric vector of the same length as the input with the probability
#' of each model relative to each other.
#' @keywords internal
.compareAIC <- function(aics) {

  probs <- vapply(aics, FUN.VALUE = c(1), FUN = function(x){
    (exp((min(aics) - x)/2))/(sum(exp((min(aics)- aics)/2)))
  })
  return(probs)

}
