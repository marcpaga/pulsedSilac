#' @rdname calculateIsotopeFraction
#' @name calculateIsotopeFraction
#' @title Calculates the incorporated isotope fraction
#'
#' @description Calculates the fraction of an isotope ratio using the following
#' formula:
#' \ifelse{html}{\out{<img src="https://latex.codecogs.com/svg.latex?\inline&space;Isotope&space;fraction_{A}&space;=&space;\frac{ratio}{ratio&space;&plus;&space;1}" title="Isotope fraction_{A} = \frac{ratio}{ratio + 1}" />}}{\deqn{Isotope fraction_{A} = \frac{ratio}{ratio + 1}}}
#'
#' The ratio should be calculated as:
#' \ifelse{html}{\out{<img src="https://latex.codecogs.com/svg.latex?\inline&space;ratio&space;=&space;isotope_{A}/isotope_{B}" title="ratio = isotope_{A}/isotope_{B}" />}}{\deqn{ratio = isotope_{A}/isotope_{B}}}
#'
#' @param x A \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or a
#' \code{SilacProteomicsExperiment} object.
#' @param ratioAssay A \code{character} with the assay name that has the ratio
#' data.
#' @param oldIsoAssay A \code{character} with the assay name that has the new
#' isotope intensity data.
#' @param newIsoAssay A \code{character} with the assay name that has the old
#' isotope intensity data.
#' @param earlyTimepoints A \code{numeric} indicating which timepoints should
#' be considered early.
#' @param lateTimepoints A \code{numeric} indicating which timepoints should
#' be considered late.
#' @param conditionCol A \code{character} indicating which column of the colData
#' data.frame indicates the different conditions.
#' @param ... Unused.
#'
#' @details If oldIsoAssay and newIsoAssay arguments are given, then the
#' earlyTimepoints and lateTiempoints arguments can be used. These can be used
#' for example if certain proteins do not have any new isotope intensity during
#' the early timepoint. Because of that, no ratio can be calculated and could
#' lead to additional missing values. If old isotope intensity is detected, then
#' a fraction of 0 for new isotope is given. Same principle applies for the
#' late timepoint but with the isotopes in reverse.
#'
#' @return a \code{SilacProteinExperiment}, \code{SilacPeptideExperiment} or
#' \code{SilacProteomicsExperiment} object with additional assays named
#' "fraction".
#'
#' @examples
#' data('wormsPE')
#' calculateIsotopeFraction(wormsPE)
#'
#' @export
setGeneric('calculateIsotopeFraction', function(x, ...){
  standardGeneric('calculateIsotopeFraction')
})

#' @rdname calculateIsotopeFraction
#' @export
setMethod('calculateIsotopeFraction', 'SilacProteinExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol) {

  ## simple fraction calculation -----------------------------------------------
  ## ratio/(ratio + 1)
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

  ## advanced fraction calculation with intensity assays -----------------------

  ## argument check ------------------------------------------------------------
  if (!oldIsoAssay %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', oldIsoAssay)
    stop(txt)
  }
  if (!newIsoAssay %in% names(assays(x))) {
    txt <- sprintf('%s not found in assay names', newIsoAssay)
    stop(txt)
  }
  if (!missing(conditionCol)) {
    metadata(x)[['conditionCol']] <- conditionCol
  }

  ## data processing -----------------------------------------------------------
  old_int <- assays(x)[[oldIsoAssay]]
  new_int <- assays(x)[[newIsoAssay]]
  total_int <- old_int + new_int


  ## loop for each condition since the early and late timepoints are relative
  ## to them
  loopCols <- .loopWrapper(x, 'conditionCol')

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
          new_fraction[missing, j] <- 1
        }
      }
    }
    assayList[[i]] <- new_fraction
  }

  fraction <- do.call('cbind', assayList)

  assays(x)[['fraction']] <- fraction
  return(x)

})

#' @rdname calculateIsotopeFraction
#' @export
setMethod('calculateIsotopeFraction', 'SilacPeptideExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol) {

  callNextMethod()

})

#' @rdname calculateIsotopeFraction
#' @export
setMethod('calculateIsotopeFraction', 'SilacProteomicsExperiment',
          function(x,
                   ratioAssay = 'ratio',
                   oldIsoAssay,
                   newIsoAssay,
                   earlyTimepoints,
                   lateTimepoints,
                   conditionCol) {


  if (any(missing(oldIsoAssay), missing(newIsoAssay))) {
    x@SilacProteinExperiment <- calculateIsotopeFraction(x@SilacProteinExperiment,
                                                    ratioAssay = ratioAssay)

    x@SilacPeptideExperiment <- calculateIsotopeFraction(x@SilacPeptideExperiment,
                                                    ratioAssay = ratioAssay)

    return(x)
  }

  x@SilacProteinExperiment <- calculateIsotopeFraction(x@SilacProteinExperiment,
                                          oldIsoAssay = oldIsoAssay,
                                          newIsoAssay = newIsoAssay,
                                          earlyTimepoints = earlyTimepoints,
                                          lateTimepoints = lateTimepoints,
                                          conditionCol = conditionCol)

  x@SilacPeptideExperiment <- calculateIsotopeFraction(x@SilacPeptideExperiment,
                                          oldIsoAssay = oldIsoAssay,
                                          newIsoAssay = newIsoAssay,
                                          earlyTimepoints = earlyTimepoints,
                                          lateTimepoints = lateTimepoints,
                                          conditionCol = conditionCol)
  return(x)

})
