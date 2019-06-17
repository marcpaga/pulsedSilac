#' Calculates the incorporated isotope fraction
#'
#' Calculates the fraction of an isotope ratio using the followin formula:
#'
#' \deqn{Isotope fraction_{A} = \frac{ratio}{ratio + 1}}
#'
#' The ratio should be calculated as:
#' \deqn{ratio = isotope_{A}/isotope_{B}}
#'
#' @param x a \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object.
#' @param assayName a \code{character} with the name of the assay that contains
#' isotope ratios.
#'
#' @return a \code{ProteinExperiment}, \code{PeptideExperiment} or
#' \code{ProteomicsExperiment} object with additional assays named "fraction".
#' @examples
#' calculateIsotopeFraction(wormsPE)
#'
#' @seealso \code{\link{ProteomicsExperiment-class}}
#' @export
setGeneric('calculateIsotopeFraction', function(x, ...){
  standardGeneric('calculateIsotopeFraction')
})

#' @param x
#' @param ratioAssay
#' @param oldIsoAssay
#' @param newIsoAssay
#' @param earlyTimepoints
#' @param lateTimepoints
#' @param conditionCol
#' @param replicateTimeCol
#' @export
setMethod('calculateIsotopeFraction', 'ProteinExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol,
                   replicateTimeCol) {

  if (any(missing(oldIsoAssay), missing(newIsoAssay))) {

    ## argument check
    if (!ratioAssay %in% names(assays(x))) {
      txt <- sprintf('%s not found in assay names', ratioAssay)
      stop(txt)
    }
    fraction_assay <- assays(x)[[ratioAssay]]/( 1 + assays(x)[[ratioAssay]] )
    assays(x)[['fraction']] <- fraction_assay

    return(x)
  }

  ## argument check
  if (!oldIsoAssay %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', oldIsoAssay)
    stop(txt)
  }
  if (!newIsoAssay %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', newIsoAssay)
    stop(txt)
  }

  old_int <- assays(x)[[oldIsoAssay]]
  new_int <- assays(x)[[newIsoAssay]]
  total_int <- old_int + new_int

  if (!missing(conditionCol)) {
    metaoptions(x)[['conditionCol']] <- conditionCol
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

  for (i in seq_along(loopCols)) {

    if (i == 1){
      assayList <- list()
    }

    ## calculate the fraction
    temp_old_int <- old_int[, loopCols[[i]]]
    temp_new_int <- new_int[, loopCols[[i]]]
    temp_total_int <- total_int[, loopCols[[i]]]

    new_fraction <- temp_new_int/temp_total_int

    ## loop over the early and late timepoints and put 0 if they have expression
    ## in one isotope and not in the other
    if (!missing(earlyTimepoints)) {
      for (j in earlyTimepoints) {

        missing <- which(is.na(temp_new_int[, j]) & !is.na(temp_old_int[, j]))
        if (length(missing) != 0) {
          new_fraction[missing, j] <- 0
        }
      }
    }

    if (!missing(lateTimepoints)) {
      for (j in lateTimepoints) {

        missing <- which(!is.na(temp_new_int[, j]) & is.na(temp_old_int[, j]))
        if (length(missing) != 0) {
          new_fraction[missing, j] <- 0
        }
      }
    }
    assayList[[i]] <- new_fraction
  }

  new_fraction <- do.call('cbind', assayList)

  assays(x)[['new_fraction']] <- new_fraction
  return(x)

})

#' @export
setMethod('calculateIsotopeFraction', 'PeptideExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol,
                   replicateTimeCol) {

  callNextMethod()

})

#' @export
setMethod('calculateIsotopeFraction', 'ProteomicsExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol,
                   replicateTimeCol) {


  if (any(missing(oldIsoAssay), missing(newIsoAssay))) {
    x@ProteinExperiment <- calculateIsotopeFraction(x@ProteinExperiment,
                                                    ratioAssay = ratioAssay)

    x@PeptideExperiment <- calculateIsotopeFraction(x@PeptideExperiment,
                                                    ratioAssay = ratioAssay)

    return(x)
  }

  x@ProteinExperiment <- calculateIsotopeFraction(x@ProteinExperiment,
                                          oldIsoAssay = oldIsoAssay,
                                          newIsoAssay = newIsoAssay,
                                          earlyTimepoints = earlyTimepoints,
                                          lateTimepoints = lateTimepoints,
                                          conditionCol = conditionCol,
                                          replicateTimeCol = replicateTimeCol)

  x@PeptideExperiment <- calculateIsotopeFraction(x@PeptideExperiment,
                                          oldIsoAssay = oldIsoAssay,
                                          newIsoAssay = newIsoAssay,
                                          earlyTimepoints = earlyTimepoints,
                                          lateTimepoints = lateTimepoints,
                                          conditionCol = conditionCol,
                                          replicateTimeCol = replicateTimeCol)
  return(x)

})
