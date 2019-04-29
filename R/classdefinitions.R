###### PROTEIN EXPERIMENT ======================================================

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.ProteinExperiment <- setClass(Class = 'ProteinExperiment',
                               contains = 'SummarizedExperiment'
)

.valid.ProteinExperiment.metaoptions<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metadata(x)))) {
    missing_names_pos <- which(!metaopions_names %in% names(metadata(x)))
    missing_names <- metaopions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metadata, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metadata(x)) %in% metaoptions_names) > 4) {

    names_pos <- which(names(metadata(x)) %in% metaoptions_names)
    names_met <- names(metadata(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metadata, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

  ## all validity is passed
  return(NULL)
}


## Wrapper for all the validity check functions
.valid.ProteinExperiment <- function(x) {

  c(.valid.ProteinExperiment.metaoptions(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteinExperiment', .valid.ProteinExperiment)

###### PEPTIDE EXPERIMENT ======================================================

#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @export
.PeptideExperiment <- setClass(Class = 'PeptideExperiment',
                               contains = 'SummarizedExperiment'
)


.valid.PeptideExperiment.metaoptions<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metadata(x)))) {
    missing_names_pos <- which(!metaopions_names %in% names(metadata(x)))
    missing_names <- metaopions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metadata, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metadata(x)) %in% metaoptions_names) > 4) {

    names_pos <- which(names(metadata(x)) %in% metaoptions_names)
    names_met <- names(metadata(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metadata, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

  ## all validity is passed
  return(NULL)
}


## Wrapper for all the validity check functions
.valid.PeptideExperiment <- function(x) {

  c(.valid.PeptideExperiment.metaoptions(x))

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
                         metadata = 'list')

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
                         'idColPep',
                         'linkedSubset',
                         'subsetMode')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metadata(x)))) {
    missing_names_pos <- which(!metaopions_names %in% names(metadata(x)))
    missing_names <- metaopions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metadata, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## duplicated metaoptions entries
  if (sum(names(metadata(x)) %in% metaoptions_names) > 8) {

    names_pos <- which(names(metadata(x)) %in% metaoptions_names)
    names_met <- names(metadata(x))[names_pos]
    dup <- names_met[which(table(names_met) > 1)]

    txt <- sprintf(
      'Duplicated metadata, the following are duplicated: %s',
      paste(dup, collapse = ' ')
    )
    return(txt)
  }

}

## Wrapper for all the validity check functions
.valid.ProteomicsExperiment <- function(x) {

  c(.valid.ProteomicsExperiment.colData(x),
    .valid.ProteomicsExperiment.metaoptions(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteomicsExperiment', .valid.ProteomicsExperiment)

