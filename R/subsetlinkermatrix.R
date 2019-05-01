#' This function subsets a linkerDf while preserving the correct links
#'
#' @param x linkerDf
#' @param protRows vector indicating which protein rows to be subsetted
#' @param protIDs vector indicating which protein IDs to be subsetted
#' @param pepRows vector indicating which peptide rows to be subsetted
#' @param pepIDs vector indicating which peptide IDs to be subsetted
#' @keywords internal
subsetLinkerDf <- function(x, protRows, protIDs, pepRows, pepIDs) {

  if (missing(x)) {
    stop('A linkerDf is needed')
  }

  if (sum(!missing(protRows),
          !missing(protIDs),
          !missing(pepRows),
          !missing(pepIDs)) > 1) {

    stop('Give only one argument')
  }

  if (!missing(protRows)) {

    ## initial subset by matcht
    i <- which(x[,'protRow'] %in% protRows)
    ## get the peptide IDs that are there
    otherPepsIDs <- unique(x[i, 'pepID'])
    ## get all the connections of those peptides too
    finalRows <- which(x[, 'pepID']  %in% otherPepsIDs)
    new.x <- x[finalRows, , drop = FALSE]

    return(new.x)
  }

  if (!missing(protIDs)) {

    ## initial subset by matcht
    i <- which(x[,'protID'] %in% protIDs)
    ## get the peptide IDs that are there
    otherPepsIDs <- unique(x[i, 'pepID'])
    ## get all the connections of those peptides too
    finalRows <- which(x[, 'pepID']  %in% otherPepsIDs)
    new.x <- x[finalRows, , drop = FALSE]

    return(new.x)
  }

  if (!missing(pepRows)) {

    ## initial subset by matcht
    i <- which(x[,'pepRow'] %in% pepRows)
    ## get the peptide IDs that are there
    otherProtIDs <- unique(x[i, 'protID'])
    ## get all the connections of those peptides too
    finalRows <- which(x[, 'protID']  %in% otherProtIDs)
    new.x <- x[finalRows, , drop = FALSE]

    return(new.x)

  }

  if (!missing(pepIDs)) {

    ## initial subset by matcht
    i <- which(x[,'pepID'] %in% pepIDs)
    ## get the peptide IDs that are there
    otherProtIDs <- unique(x[i, 'protID'])
    ## get all the connections of those peptides too
    finalRows <- which(x[, 'protID']  %in% otherProtIDs)
    new.x <- x[finalRows, , drop = FALSE]

    return(new.x)

  }

}
