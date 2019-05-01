###### Logical checks ==========================================================

setMethod('hasRowData', 'ProteinExperiment', function(x){

  rd <- rowData(x)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

})

setMethod('hasRowData', 'PeptideExperiment', function(x){

  rd <- rowData(x)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

})

setMethod('hasRowData', 'ProteomicsExperiment', function(x){

  outVec <- c(NA, NA)
  rd <- rowData(x@ProteinExperiment)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    outVec[1] <- TRUE
  } else {
    outVec[1] <- FALSE
  }

  rd <- rowData(x@PeptideExperiment)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    outVec[2] <- TRUE
  } else {
    outVec[2] <- FALSE
  }

  return(outVec)
})




###### Specific metaoptions ====================================================

## retrieve metaoptions =====

hasMetaoption <- function(x, option) {

  op <- metaoptions(x)[[option]]
  if (is.na(op)) {
    txt <- sprintf('Not defined in metadata: %s.', option)
    stop(txt)
  } else {
    return(TRUE)
  }

}

metaoptionInColData <- function(x, option) {

  if (hasMetaoption(x, option)) {
    op <- metaoptions(x)[[option]]

    if (op %in% names(colData(x))) {
      return(TRUE)
    } else {
      txt <- sprintf('Column not found in colData: %s.', op)
      stop(txt)
    }
  }

}

giveMetaoption <- function(x, option) {

  if (metaoptionInColData(x, option)) {
    op <- metaoptions(x)[[option]]
    return(op)
  }

}

## merge metaoptions lists ====

.mergeMetaoptions <- function(x, y) {

  finalLength <- max(length(x), length(y))
  finalNames <- unique(c(names(x), names(y)))
  new_metaoptions <- list()

  for (i in seq_len(finalLength)) {

    new_metaoptions[[i]] <- .compareOptions(x[[finalNames[i]]],
                                            y[[finalNames[i]]],
                                            finalNames[i])

  }

  names(new_metaoptions) <- finalNames
  return(new_metaoptions)
}

.compareOptions <- function(x, y, option) {

  ## one of the elements is not in the other list
  if (is.null(x)) {
    return(y)
  }
  if (is.null(y)) {
    return(x)
  }

  ## not defined in any of the lists
  if (is.na(x) & is.na(y)) {
    return(NA)
  }

  ## defined in one but not in the other
  if (is.na(x) & !is.na(y)) {
    return(y)
  }
  if (!is.na(x) & is.na(y)) {
    return(x)
  }

  ## defined in both
  if (x == y) {
    return(x)
  } else if (x != y) {
    txt <- sprintf(
      paste('The %s metaoption was defined differently: %s and %s in the two',
            'experiments. The first one was used.'),
      option, x, y
    )
    warning(txt)
    return(x)
  }
}

###### Specific linkerDf

## check integrity of the linkerDf

#' @keywords internal
checkLinkerDf <- function(m){

  # peptides have links to proteins not in the data
  missingLinksPepToProt <- sum(is.na(m[,3]))
  # proteins have links to peptides not in the data
  missingLinksProtToPep <- sum(is.na(m[,4]))

  # there are no missing links
  if (sum(missingLinksPepToProt, missingLinksProtToPep) == 0) {
    return(m)
  }

  # missing links are removed
  if (missingLinksPepToProt > 0) {

    m <- m[-which(is.na(m[,3])),]
    message(paste0(missingLinksPepToProt,
                   ' links between given peptides are not matched to proteins'))

  }

  # check again if there are still missing links since they could have already
  # been removed in the previous step
  missingLinksProtToPep <- sum(is.na(m[,4]))
  # missing links are removed
  if (missingLinksProtToPep > 0) {

    m <- m[-which(is.na(m[,4])),]
    message(paste0(missingLinksProtToPep,
                   ' links between given proteins are not matched to peptides'))

  }

  return(m)

}
