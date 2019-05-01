###### PROTEIN EXPERIMENT ======================================================

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.ProteinExperiment <- setClass(Class = 'ProteinExperiment',
                               slots = representation(metaoptions = 'list'),
                               contains = 'SummarizedExperiment'
)

.valid.ProteinExperiment.metaoptions <- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metaoptions(x)))) {
    missing_names_pos <- which(!metaoptions_names %in% names(metaoptions))
    missing_names <- metaoptions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metaoptions, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metaoptions(x)) %in% metaoptions_names) > 4) {

    names_pos <- which(names(metaoptions(x)) %in% metaoptions_names)
    names_met <- names(metaoptions(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metaoptions, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

  ## metaoptions have to be NA, numeric or character
  ## checking if the column is at the colData/rowData is done on function call
  ## not on validity because it could cause problems when changing those
  ## DataFrames

  for (mopt in metaoptions_names) {
    val <- metaoptions(x)[[mopt]]
    if (any(is.na(val), is.character(val), is.numeric(val))) {
      next
    } else {
      txt <- sprintf('%s must be character, numeric or NA', mopt)
      return(txt)
    }
  }

  ## all validity is passed
  return(NULL)
}

.valid.ProteinExperiment.rowData <- function(x) {

  if (is.null(x@elementMetadata)) {
    return(NULL)
  }

  rd_rows <- nrow(x@elementMetadata)
  if (rd_rows != nrow(x)) {
    txt <- sprintf(
      paste('The number of rows of rowData (%d) does not match with the',
            'number of rows of the ProteinExperiment (%d)'),
      rd_rows, nrow(x))
    return(txt)
  }

  return(NULL)

}


## Wrapper for all the validity check functions
.valid.ProteinExperiment <- function(x) {

  c(.valid.ProteinExperiment.metaoptions(x),
    .valid.ProteinExperiment.rowData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteinExperiment', .valid.ProteinExperiment)

###### PEPTIDE EXPERIMENT ======================================================

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @export
.PeptideExperiment <- setClass(Class = 'PeptideExperiment',
                               slots = representation(metaoptions = 'list'),
                               contains = 'SummarizedExperiment'
)


.valid.PeptideExperiment.metaoptions<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol',
                         'proteinCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metaoptions(x)))) {
    missing_names_pos <- which(!metaoptions_names %in% names(metaoptions(x)))
    missing_names <- metaoptions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metaoptions, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metaoptions(x)) %in% metaoptions_names) > 5) {

    names_pos <- which(names(metaoptions(x)) %in% metaoptions_names)
    names_met <- names(metaoptions(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metaoptions, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

  ## metaoptions have to be NA, numeric or character
  ## checking if the column is at the colData/rowData is done on function call
  ## not on validity because it could cause problems when changing those
  ## DataFrames

  for (mopt in metaoptions_names) {
    val <- metaoptions(x)[[mopt]]
    if (any(is.na(val), is.character(val), is.numeric(val))) {
      next
    } else {
      txt <- sprintf('%s must be character, numeric or NA', mopt)
      return(txt)
    }
  }

  ## all validity is passed
  return(NULL)
}

.valid.PeptideExperiment.rowData <- function(x) {

  if (is.null(x@elementMetadata)) {
    return(NULL)
  }

  rd_rows <- nrow(x@elementMetadata)
  if (rd_rows != nrow(x)) {
    txt <- sprintf(
      paste('The number of rows of rowData (%d) does not match with the',
            'number of rows of the PeptideExperiment (%d)'),
      rd_rows, nrow(x))
    return(txt)
  }

  return(NULL)

}

## Wrapper for all the validity check functions
.valid.PeptideExperiment <- function(x) {

  c(.valid.PeptideExperiment.metaoptions(x),
    .valid.PeptideExperiment.rowData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('PeptideExperiment', .valid.PeptideExperiment)

###### PROTEOMICS EXPERIMENT ===================================================

#' @title  ProteomicsExperiment class
#'
#' @export
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import methods
.ProteomicsExperiment <- setClass('ProteomicsExperiment',

                                  slots = representation(ProteinExperiment = 'ProteinExperiment',
                                                         PeptideExperiment = 'PeptideExperiment',
                                                         colData = 'DataFrame',
                                                         linkerDf = 'data.frame',
                                                         metadata = 'list',
                                                         metaoptions = 'list')

)


.valid.ProteomicsExperiment.colData <- function(x) {

  ncol_prot <- ncol(x@ProteinExperiment)
  ncol_pept <- ncol(x@PeptideExperiment)
  if (ncol_prot != ncol_pept) {
    txt <- sprintf(
      paste('\n Number of columns of ProteinExperiment (%d) must be equal to',
            'the number of columns of PeptideExperiment (%d)'),
      ncol_prot, ncol_pept
    )
    return(txt)
  } else {
    return(NULL)
  }

}

.valid.ProteomicsExperiment.metaoptions <- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol',
                         'idColProt',
                         'idColPept',
                         'linkedSubset',
                         'subsetMode',
                         'proteinCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metaoptions(x)))) {
    missing_names_pos <- which(!metaoptions_names %in% names(metaoptions(x)))
    missing_names <- metaoptions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metaoptions, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metaoptions(x)) %in% metaoptions_names) > 9) {

    names_pos <- which(names(metaoptions(x)) %in% metaoptions_names)
    names_met <- names(metaoptions(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metaoptions, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

  ## metaoptions have to be NA, numeric or character
  ## checking if the column is at the colData/rowData is done on function call
  ## not on validity because it could cause problems when changing those
  ## DataFrames

  for (mopt in metaoptions_names) {
    val <- metaoptions(x)[[mopt]]

    ## hardcode two special metaoptions
    if (mopt == 'subsetMode') {
      if (val %in% c('protein', 'peptide')) {
        next
      } else {
        txt <- 'subsetMode must be a character: "protein" or "peptide"'
        return(txt)
      }
    }

    if (mopt == 'linkedSubset') {
      if (val %in% c(TRUE, FALSE)) {
        next
      } else {
        txt <- 'linkedSubset must be a logical: TRUE or FALSE'
        return(txt)
      }
    }

    ## the rest is all either numeric, character or NA
    if (any(is.na(val), is.character(val), is.numeric(val))) {
      next
    } else {
      txt <- sprintf('%s must be character, numeric or NA', mopt)
      return(txt)
    }
  }

}

.valid.ProteomicsExperiment.linkerDf <- function(x){

  linkM <- x@linkerDf

  # no linker matrix, no need to check
  if (nrow(linkM) == 0) {
    return(NULL)
  }

  if (!is.data.frame(linkM)) {
    return('The linkerDf should be a data.frame')
  }

  if (ncol(linkM) != 4) {
    return('The linkerDf should have 4 columns')
  }

  if (x@metaoptions[['linkedSubset']]) {

    if (length(unique(linkM[, 3])) < nrow(rowDataProt(x))) {
      return('Not all proteins have a link')
    }

    if (length(unique(linkM[, 4])) < nrow(rowDataPept(x))) {
      return('Not all peptides have a link')
    }
  }

  return(NULL)

}

## Wrapper for all the validity check functions
.valid.ProteomicsExperiment <- function(x) {

  c(.valid.ProteomicsExperiment.colData(x),
    .valid.ProteomicsExperiment.metaoptions(x),
    .valid.ProteomicsExperiment.linkerDf(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteomicsExperiment', .valid.ProteomicsExperiment)

