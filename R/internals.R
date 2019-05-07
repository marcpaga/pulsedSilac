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

#' @keywords internal
hasAssayNames <- function(x) {

  if (is.null(assayNames(x))) {
    return(FALSE)
  } else {
    return(TRUE)
  }

}


###### Specific metaoptions ====================================================

## retrieve metaoptions =====

#' @keywords internal
hasMetaoption <- function(x, option) {

  op <- metaoptions(x)[[option]]
  if (is.na(op)) {
    txt <- sprintf('Not defined in metaoptions: %s.', option)
    stop(txt)
  } else {
    return(TRUE)
  }

}

#' @keywords internal
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

#' @keywords internal
metaoptionInRowData <- function(x, option) {

  if (hasMetaoption(x, option)) {
    op <- metaoptions(x)[[option]]

    if (option == 'idColProt' & op %in% colnames(rowDataProt(x))) {
      return(TRUE)
    } else if (option == 'idColPept' & op %in% colnames(rowDataPept(x))) {
      return(TRUE)
    } else {
      txt <- sprintf('Column not found in colData: %s.', op)
      stop(txt)
    }
  }

}

#' @keywords internal
giveMetaoption <- function(x, option) {

  ## hardcode for ProteomicsExperiment unique options
  if (option == 'subsetMode') {
    return(metaoptions(x)[[option]])
  }

  if (option == 'linkedSubset') {
    return(metaoptions(x)[[option]])
  }

  if (option %in% c('idColProt', 'idColPept')) {
    if (metaoptionInRowData(x, option)) {
      op <- metaoptions(x)[[option]]
      return(op)
    }
  }

  if (metaoptionInColData(x, option)) {
    op <- metaoptions(x)[[option]]
    return(op)
  }

}

## merge metaoptions lists ====

#' @keywords internal
mergeMetaoptions <- function(x, y) {

  finalLength <- max(length(x), length(y))
  finalNames <- unique(c(names(x), names(y)))
  new_metaoptions <- list()

  for (i in seq_len(finalLength)) {

    new_metaoptions[[i]] <- compareOptions(x[[finalNames[i]]],
                                           y[[finalNames[i]]],
                                           finalNames[i])

  }

  names(new_metaoptions) <- finalNames
  return(new_metaoptions)
}

#' @keywords internal
compareOptions <- function(x, y, option) {

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

## synchronize metaoptions =====

#' This function sends the values from the metaoptions of ProteomicsExperiment
#' to the ProteinExperiment and PeptideExperiment
#' @keywords internal
synchronizeMetaoptions <- function(x) {

  for (meta in names(metaoptions(x))) {

    value <- metaoptions(x)[[meta]]
    if (meta %in% names(metaoptions(x@ProteinExperiment))) {
      metaoptions(x@ProteinExperiment)[[meta]] <- value
    }

    if (meta %in% names(metaoptions(x@PeptideExperiment))) {
      metaoptions(x@PeptideExperiment)[[meta]] <- value
    }
  }
  validObject(x)
  return(x)

}

###### Specific linkerDf -----

## check integrity of the linkerDf----

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

## rbind two linker data frames

#' @keywords internal
rbindLinkerDf <- function(x, y) {


  ## need to remove NA, they can be a results of a subset function
  if (sum(is.na(x[,3])) > 0) {
    x <- x[-which(is.na(x[,3])),]
  }
  if (sum(is.na(y[,3])) > 0) {
    y <- y[-which(is.na(y[,3])),]
  }

  proteinIntersect <- intersect(x[which(!is.na(x[, 3])), 1],
                                y[which(!is.na(y[, 3])), 1])
  peptideIntersect <- intersect(x[which(!is.na(x[, 4])), 2],
                                y[which(!is.na(y[, 4])), 2])

  temp.pr <- subset(x, !is.na(x[, 3]))
  max.pr <- length(unique(temp.pr[, 3]))
  temp.pe <- subset(x, !is.na(x[, 4]))
  max.pe <- length(unique(temp.pe[, 4]))


  ## different protein IDs and different peptide IDs
  if (length(proteinIntersect) == 0 & length(peptideIntersect) == 0) {
    new.y <- y
    new.y[, 3] <- max.pr + y[, 3]
    new.y[, 4] <- max.pe + y[, 4]

    new.lm <- rbind(x, new.y)
    rownames(new.lm) <- seq_len(nrow(new.lm))
    return(new.lm)
  }

  ## protein overlap, but different peptide IDs
  if (length(proteinIntersect) != 0 & length(peptideIntersect) == 0) {

    new.y <- y
    for (row in seq_len(nrow(new.y))) {
      e <- new.y[row, 1]
      if (e %in% proteinIntersect){
        new.row <- unique(x[which(x[, 1] == e), 3])
        if (!is.na(new.row)) {
          new.y[which(y[, 1] == e), 3] <- new.row
        } else {
          new.y[row, 3] <- y[row, 3] + max.pr - length(proteinIntersect)
        }
        new.y[row, 4] <- y[row, 4] + max.pe
      } else {
        new.y[row, 3] <- y[row, 3] + max.pr - length(proteinIntersect)
        new.y[row, 4] <- y[row, 4] + max.pe
      }

    }

    new.lm <- rbind(x, new.y)
    rownames(new.lm) <- seq_len(nrow(new.lm))
    return(new.lm)
  }

  ## different protein IDs, but overlap in peptide IDs
  if (length(proteinIntersect) == 0 & length(peptideIntersect) != 0) {

    new.y <- y
    for (row in seq_len(nrow(new.y))) {

      e <- new.y[row, 2]
      if (e %in% peptideIntersect){
        new.row <- unique(x[which(x[, 2] == e), 4])
        if (!is.na(new.row)) {
          new.y[which(y[, 2] == e), 4] <- new.row
        } else {
          new.y[row, 4] <- y[row, 4] + max.pe - length(peptideIntersect)
        }

        new.y[row, 3] <- y[row, 3] + max.pr
      } else {
        new.y[row, 3] <- y[row, 3] + max.pr
        new.y[row, 4] <- y[row, 4] + max.pe - length(peptideIntersect)
      }

    }

    new.lm <- rbind(x, new.y)
    rownames(new.lm) <- seq_len(nrow(new.lm))
    return(new.lm)
  }

  ## protein ID overlap and peptide ID overlap
  if (length(proteinIntersect) != 0 & length(peptideIntersect) != 0) {

    new.y <- y
    for (row in seq_len(nrow(new.y))) {

      e <- new.y[row, 1]
      if (e %in% proteinIntersect){
        new.row <- unique(x[which(x[, 1] == e), 3])
        if (!is.na(new.row)) {
          new.y[row, 3] <- new.row
        } else {
          new.y[row, 3] <- y[row, 3] + max.pr - length(proteinIntersect)
        }
      } else {
        new.y[row, 3] <- y[row, 3] + max.pr - length(proteinIntersect)
      }

      e <- new.y[row, 2]
      if (e %in% peptideIntersect){
        new.row <- unique(x[which(x[, 2] == e), 4])
        if (!is.na(new.row)) {
          new.y[row, 4] <- new.row
        } else {
          new.y[row, 4] <- y[row, 4] + max.pe - length(peptideIntersect)
        }
      } else {
        new.y[row, 4] <- y[row, 4] + max.pe - length(peptideIntersect)
      }

    }

    new.lm <- rbind(x, new.y)
    if (sum(duplicated(new.lm)) > 0) {
      new.lm <- new.lm[-which(duplicated(new.lm)),]
    }

    rownames(new.lm) <- seq_len(nrow(new.lm))
    return(new.lm)
  }

}


###### Experiment looping ======================================================

## wrapper for the following 3 functions

#' @keywords internal
experimentLoopWrapper <- function(x, type) {

  proCD <- processColData(x, type = type)
  expCols <- getExperimentCols(proCD)
  loopCols <- getLoopCols(expCols, type = type)
  return(loopCols)

}


## The following are 3 functions called in the following order:
##    1. processColData
##    2. getExperimentCols
##    3. getLoopCols
## This functions process the colData from a SummarizedExperiment-like object
## together with the metaoptions and return a list of columns to loop over.
## The columns are different depending on the analysis that is being done.

### processColData ====
processColData <- function(x, type) {

  df <- colData(x)

  ## looping that goes over each experiment condition and time replicates
  ## for this type we need 'conditionCol' and 'replicateTimeCol'
  if (type == 'cond.timerep') {

    conditionCol <- giveMetaoption(x = x, option = 'conditionCol')
    replicateTimeCol <- giveMetaoption(x = x, option = 'replicateTimeCol')

    i <- which(colnames(colData(x)) == conditionCol)
    colnames(df)[i] <- 'conditionCol'
    i <- which(colnames(colData(x)) == replicateTimeCol)
    colnames(df)[i] <- 'replicateTimeCol'

  } else if (type == 'cond') {

    conditionCol <- giveMetaoption(x = x, option = 'conditionCol')

    i <- which(colnames(colData(x)) == conditionCol)
    colnames(df)[i] <- 'conditionCol'

  }

  df <- data.frame(lapply(df, as.factor))

  i <- which(colnames(df) %in% c('conditionCol', 'timeCol',
                                 'replicateIntCol', 'replicateTimeCol'))
  df <- df[, i, drop = FALSE]

  return(df)

}

### getExperimentCols ====
getExperimentCols <- function(x) {


  experimentCols <- as.list(x)
  if (!'conditionCol' %in% names(experimentCols)) {
    experimentCols$conditionCol <- factor(NULL)
  }
  if (!'timeCol' %in% names(experimentCols)) {
    experimentCols$timeCol <- factor(NULL)
  }
  if (!'replicateTimeCol' %in% names(experimentCols)) {
    experimentCols$replicateTimeCol <- factor(NULL)
  }
  if (!'replicateIntCol' %in% names(experimentCols)) {
    experimentCols$replicateIntCol <- factor(NULL)
  }

  return(experimentCols)

}

### getLoopCols ====
getLoopCols <- function(x, type) {

  listnames <- c('conditionCol', 'timeCol',
                 'replicateIntCol', 'replicateTimeCol')
  if (!all(listnames %in% names(x))) {
    stop('Input must be a list with the following names: conditionalCol, ',
         'timeCol, replicateIntCol and replicateTimeCol')
  }

  if (type == 'cond.timerep') {

    samples <- factor(paste0(x[['conditionCol']], x[['replicateTimeCol']]))
    cols <- lapply(levels(samples), function(x) which(samples == x))
    return(cols)

  } else if (type == 'cond') {

    samples <- factor(x[['conditionCol']])
    cols <- lapply(levels(samples), function(x) which(samples == x))
    return(cols)
  }

  stop('Please provide a type to group the columns')

}

