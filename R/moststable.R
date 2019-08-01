#' @rdname mostStable
#' @name mostStable
#' @title Most stable proteins/peptides
#'
#' @description Finds which are the most stable proteins/peptides across
#' the entire experiment. These proteins/peptides can be used to estimate
#' the cell growth of each condition.
#'
#' @details Proteins/peptides that are not found in all timepoints and in all
#' conditions are not considered. The stability is based on ranking
#' heavy label incorporation for each timepoint; therefore, lower values
#' are correlated to higher stability.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param assayName Name of the assay to use.
#' @param n A \code{numeric} indicating how many proteins should be returned.
#' @param mode A \code{character} indicating which level of data to use,
#' either "protein" or "peptide". Only relevant for ProteomicsExperiment
#' inputs.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param replicateTimeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different time replicates.
#' @param ... Unused.
#'
#' @return A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object with the n most stable proteins/peptides.
#'
#' @examples
#' mostStable(mefPE, assayName = 'fraction', n = 50)
#' @export
setGeneric('mostStable', function(x, ...){
  standardGeneric('mostStable')
})


#' @rdname mostStable
#' @export
setMethod('mostStable',
          'ProteinExperiment',
          function(x,
                   assayName,
                   n,
                   conditionCol,
                   replicateTimeCol) {


  ## argument check
  if (!assayName %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', assayName)
    stop(txt)
  }

  mat <- assays(x)[[assayName]]

  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
  }
  if (!missing(replicateTimeCol)) {
    metaoptions(x)[['replicateTimeCol']] <- replicateTimeCol
  }

  ## which columns belong to which experiment
  loopList <- tryCatch(
    {
      loop <- experimentLoopWrapper(x, 'cond.timerep')
      loop
    },
    error = function(c){
      tryCatch(
        {
          loop <- experimentLoopWrapper(x, 'cond')
          loop
        },
        error = function(c){
          list(seq_len(ncol(x)), NA)
        }
      )
    }
  )

  rankRes <- matrix(data = NA,
                    ncol = length(loopList),
                    nrow = nrow(x))

  for (j in seq_along(loopList)) {

    tempMat <- assays(x)[[assayName]][, loopList[[j]]]
    rankMat <- apply(tempMat, 2, rank, na.last = 'keep')
    resMat <- rank(apply(rankMat, 1, sum, na.rm = FALSE), na.last = 'keep')
    rankRes[, j] <- resMat
  }

  cutoff <- sort(apply(rankRes, 1, mean))[n]
  stablePE <- x[which(apply(rankRes, 1, mean) <= cutoff), ]

  return(stablePE)

})

#' @rdname mostStable
#' @export
setMethod('mostStable',
          'PeptideExperiment',
          function(x,
                   assayName,
                   n,
                   conditionCol,
                   replicateTimeCol) {

  callNextMethod()

})

#' @rdname mostStable
#' @export
setMethod('mostStable',
          'ProteomicsExperiment',
          function(x,
                   assayName,
                   n,
                   mode,
                   conditionCol,
                   replicateTimeCol) {

  if (mode == 'protein') {
    mostStable(x = x@ProteinExperiment,
               assayName = assayName,
               n = n,
               conditionCol = conditionCol,
               replicateTimeCol = replicateTimeCol)
  } else if (mode == 'peptide') {
    mostStable(x = x@PeptideExperiment,
               assayName = assayName,
               n = n,
               conditionCol = conditionCol,
               replicateTimeCol = replicateTimeCol)
  }

})
