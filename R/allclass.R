###### PROTEIN EXPERIMENT ======================================================

#' @rdname SilacProteinExperiment-class
#' @name SilacProteinExperiment-class
#' @title  SilacProteinExperiment class
#'
#' @description S4 class that extends the \code{\link{SummarizedExperiment}}
#' class. This class is designed for proteomics data, more especifically
#' protein level data. The \code{metadata} slot comes already initialized with
#' the metaoptions (see details).
#'
#' @details The \code{SilacProteinExperiment} class has been designed to store
#' protein level data and to be used in the functions provided in this package
#' for pulsed SILAC data analysis; in combination with the other two classes
#' from the package: the \code{\link{SilacPeptideExperiment}} and
#' \code{\link{SilacProteomicsExperiment}} classes.
#'
#' SilacProteinExperiment metaoptions are stored in the \code{metadata} slot
#' This contains a \code{list} with some parameters that are automatically
#' initialized by the constructor. Some parameters are mandatory for certain
#' functions or operations. The user can add or remove items at their
#' discretion. These parameters are meant to help automate certain pipeline or
#' data analysis steps. These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'}
#'
#' @section Constructor:
#' See \link{SilacProteinExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{SilacProteinPeptideExperiment-accessors} for details.
#'
#' @seealso \code{\link{SilacProteinExperiment-constructor}},
#'          \code{\link{SilacProteinPeptideExperiment-accessors}},
#'          \code{\link{SummarizedExperiment}}
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @export
.SilacProteinExperiment <- setClass(Class = 'SilacProteinExperiment',
                               contains = 'SummarizedExperiment'
)

.valid.SilacProteinExperiment.metaoptions <- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol')

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

  ## all validity is passed
  return(NULL)
}

.valid.SilacProteinExperiment.rowData <- function(x) {

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

.valid.SilacProteinExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}


## Wrapper for all the validity check functions
.valid.SilacProteinExperiment <- function(x) {

  c(.valid.SilacProteinExperiment.metaoptions(x),
    .valid.SilacProteinExperiment.rowData(x),
    .valid.SilacProteinExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('SilacProteinExperiment', .valid.SilacProteinExperiment)

###### PEPTIDE EXPERIMENT ======================================================

#' @rdname SilacPeptideExperiment-class
#' @name SilacPeptideExperiment-class
#' @title  SilacPeptideExperiment class
#'
#' @description S4 class that extends the \code{\link{SummarizedExperiment}}
#' class. This class is designed for proteomics data, more especifically
#' peptide level data. The \code{metadata} slot comes already initialized with
#' the metaoptions (see details).
#'
#' @details The \code{SilacPeptideExperiment} class has been designed to store
#' peptide level data and to be used in the functions provided in this package
#' for pulsed SILAC data analysis; in combination with the other two classes
#' from the package: the \code{\link{SilacProteinExperiment}} and
#' \code{\link{SilacProteomicsExperiment}} classes.
#'
#' ProteinExperiment metaoptions are stored in the \code{metadata} slot
#' This contains a \code{list} with some parameters that are automatically
#' initialized by the constructor. Some parameters are mandatory for certain
#' functions or operations. The user can add or remove items at their
#' discretion. These parameters are meant to help automate certain pipeline or
#' data analysis steps. These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'   \item{proteinCol}{\code{character} indicating the column name of
#'   \code{rowData(x)} that defines to which protein a peptide is assigned.}
#'}
#'
#' @section Constructor:
#' See \link{SilacPeptideExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{SilacProteinPeptideExperiment-accessors} for details.
#'
#' @seealso \code{\link{SilacPeptideExperiment-constructor}},
#'          \code{\link{SilacProteinPeptideExperiment-accessors}},
#'          \code{\link{SummarizedExperiment}}
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.SilacPeptideExperiment <- setClass(Class = 'SilacPeptideExperiment',
                               contains = 'SilacProteinExperiment'
)


.valid.SilacPeptideExperiment.metaoptions<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'proteinCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metadata(x)))) {
    missing_names_pos <- which(!metaoptions_names %in% names(metaoptions(x)))
    missing_names <- metaoptions_names[missing_names_pos]

    txt <- sprintf(
      'Incomplete metaoptions, the following are missing: %s',
      paste(missing_names, collapse = ' ')
    )
    return(txt)
  }

  ## all validity is passed
  return(NULL)
}

.valid.SilacPeptideExperiment.rowData <- function(x) {

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

.valid.SilacPeptideExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}

## Wrapper for all the validity check functions
.valid.SilacPeptideExperiment <- function(x) {

  c(.valid.SilacPeptideExperiment.metaoptions(x),
    .valid.SilacPeptideExperiment.rowData(x),
    .valid.SilacPeptideExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('SilacPeptideExperiment', .valid.SilacPeptideExperiment)

###### PROTEOMICS EXPERIMENT ===================================================

#' @rdname SilacProteomicsExperiment-class
#' @name SilacProteomicsExperiment-class
#' @title  SilacProteomicsExperiment class
#'
#' @description S4 class that contains a \code{SilacProteinExperiment} object
#' and a \code{SilacPeptideExperiment} object. The two objects are linked by a
#' \code{data.frame} (linkerDf). This class can be used to manage both protein
#' and peptide data at the same time.
#'
#' @slot SilacProteinExperiment Contains a \code{SilacProteinExperiment} object.
#' @slot SilacPeptideExperiment Contains a \code{SilacPeptideExperiment} object.
#' @slot colData Contains a \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @slot linkerDf Contains a \code{data.frame} that has been created with
#' \code{\link{buildLinkerDf}}. It contains the relationships between proteins
#' and peptides.
#' @slot metadata Contains a \code{list} to store any kind of experiment-wide
#' data and the metaoptions.

#' @details The \code{SilacProteomicsExperiment} object is just a
#' \code{SilacProteinExperiment} object and a \code{SilacPeptideExperiment}
#' object together.
#'
#' The rows of the \code{SilacProteinExperiment} object represent proteins.
#' The rows of the \code{SilacPeptideExperiment} object represent peptides.
#'
#' The columns of the \code{SilacProteomicsExperiment} object represent samples.
#' Samples are shared at both protein and peptide levels.
#'
#' Experiment-wide information can be stored in the \code{metadata} slot, which
#' is accessed with the \code{metadata} function. This contains a \code{list}
#' object in which each item is left to the discretion of the user. Some
#' possible examples could be: data of the experiment, author, machine used,
#' etc.
#'
#' SilacProteomicsExperiment options are stored in the \code{metadata} slot.
#' This contains a \code{list} with some parameters that are automatically
#' initialized by the constructor. Some parameters are mandatory for certain
#' functions or operations. The user can add or remove items at their
#' discretion. These parameters are meant to help automate certain pipeline or
#' data analysis steps. These metaoptions are:
#' These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'   \item{idColProt}{A \code{character} indicating which column from the
#'   rowData (protein) should be used as ids.}
#'   \item{idColPept}{A \code{character} indicating which column from the
#'   rowData (peptides) should be used as ids.}
#'   \item{linkedSubset}{A \code{logical} if subsetting should be linked between
#'   proteins and peptide.}
#'   \item{subsetMode}{A \code{character}, either 'protein' or 'peptide'
#'   indicating which level should be used first when subsetting.}
#'}
#'
#' @section Constructor:
#' See \link{SilacProteomicsExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{SilacProteomicsExperiment-accessors} for details.
#'
#' @seealso \code{\link{SilacProteomicsExperiment-constructor}},
#'          \code{\link{SilacProteomicsExperiment-accessors}},
#'          \code{\link{SilacProteinExperiment}},
#'          \code{\link{SilacPeptideExperiment}}

#' @export
.SilacProteomicsExperiment <- setClass('SilacProteomicsExperiment',

  slots = representation(SilacProteinExperiment = 'SilacProteinExperiment',
                         SilacPeptideExperiment = 'SilacPeptideExperiment',
                         colData = 'DataFrame',
                         linkerDf = 'data.frame',
                         metadata = 'list')

)


.valid.SilacProteomicsExperiment.ncol <- function(x) {

  ncol_prot <- ncol(x@SilacProteinExperiment)
  ncol_pept <- ncol(x@SilacPeptideExperiment)
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

.valid.SilacProteomicsExperiment.metaoptions <- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
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
  }

}

.valid.SilacProteomicsExperiment.linkerDf <- function(x){

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

  if (metaoptions(x)[['linkedSubset']]) {

    if (length(unique(linkM[, 3])) < nrow(rowDataProt(x))) {
      return('Not all proteins have a link')
    }

    if (length(unique(linkM[, 4])) < nrow(rowDataPept(x))) {
      return('Not all peptides have a link')
    }
  }

  return(NULL)

}

.valid.SilacProteomicsExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}

## Wrapper for all the validity check functions
.valid.SilacProteomicsExperiment <- function(x) {

  c(.valid.SilacProteomicsExperiment.ncol(x),
    .valid.SilacProteomicsExperiment.metaoptions(x),
    .valid.SilacProteomicsExperiment.linkerDf(x),
    .valid.SilacProteomicsExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('SilacProteomicsExperiment', .valid.SilacProteomicsExperiment)

