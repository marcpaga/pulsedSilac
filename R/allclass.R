###### PROTEIN EXPERIMENT ======================================================

#' @rdname ProteinExperiment-class
#' @name ProteinExperiment-class
#' @title  ProteinExperiment class
#'
#' @description S4 class that extends the \code{\link{SummarizedExperiment}}
#' class. This class is designed for proteomics data, more especifically
#' protein level data. The \code{metadata} slot comes already initialized with
#' the metaoptions (see details).
#'
#' @details The \code{ProteinExperiment} class has been designed to store
#' protein level data and to be used in the functions provided in this package
#' for pulsed SILAC data analysis; in combination with the other two classes
#' from the package: the \code{\link{PeptideExperiment}} and
#' \code{\link{ProteomicsExperiment}} classes.
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
#'}
#'
#' @section Constructor:
#' See \link{ProteinExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{ProteinPeptideExperiment-accessors} for details.
#'
#' @seealso \code{\link{ProteinExperiment-constructor}},
#'          \code{\link{ProteinPeptideExperiment-accessors}},
#'          \code{\link{SummarizedExperiment}}
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @export
.ProteinExperiment <- setClass(Class = 'ProteinExperiment',
                               contains = 'SummarizedExperiment'
)

.valid.ProteinExperiment.metaoptions <- function(x) {

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

.valid.ProteinExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}


## Wrapper for all the validity check functions
.valid.ProteinExperiment <- function(x) {

  c(.valid.ProteinExperiment.metaoptions(x),
    .valid.ProteinExperiment.rowData(x),
    .valid.ProteinExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteinExperiment', .valid.ProteinExperiment)

###### PEPTIDE EXPERIMENT ======================================================

#' @rdname PeptideExperiment-class
#' @name PeptideExperiment-class
#' @title  PeptideExperiment class
#'
#' @description S4 class that extends the \code{\link{SummarizedExperiment}}
#' class. This class is designed for proteomics data, more especifically
#' peptide level data. The \code{metadata} slot comes already initialized with
#' the metaoptions (see details).
#'
#' @details The \code{PeptideExperiment} class has been designed to store
#' peptide level data and to be used in the functions provided in this package
#' for pulsed SILAC data analysis; in combination with the other two classes
#' from the package: the \code{\link{ProteinExperiment}} and
#' \code{\link{ProteomicsExperiment}} classes.
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
#' See \link{PeptideExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{ProteinPeptideExperiment-accessors} for details.
#'
#' @seealso \code{\link{PeptideExperiment-constructor}},
#'          \code{\link{ProteinPeptideExperiment-accessors}},
#'          \code{\link{SummarizedExperiment}}
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.PeptideExperiment <- setClass(Class = 'PeptideExperiment',
                               contains = 'ProteinExperiment'
)


.valid.PeptideExperiment.metaoptions<- function(x) {

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

.valid.PeptideExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}

## Wrapper for all the validity check functions
.valid.PeptideExperiment <- function(x) {

  c(.valid.PeptideExperiment.metaoptions(x),
    .valid.PeptideExperiment.rowData(x),
    .valid.PeptideExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('PeptideExperiment', .valid.PeptideExperiment)

###### PROTEOMICS EXPERIMENT ===================================================

#' @rdname ProteomicsExperiment-class
#' @name ProteomicsExperiment-class
#' @title  ProteomicsExperiment class
#'
#' @description S4 class that contains a \code{ProteinExperiment} object and
#' a \code{PeptideExperiment} object. The two objects are linked by a
#' \code{data.frame} (linkerDf). This class can be used to manage both protein
#' and peptide data at the same time.
#'
#' @slot ProteinExperiment Contains \code{ProteinExperiment} object.
#' @slot PeptideExperiment Contains \code{PeptideExperiment} object.
#' @slot colData Contains a \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @slot linkerDf Contains a \code{data.frame} that has been created with
#' \code{\link{buildLinkerDf}}. It contains the relationships between proteins
#' and peptides.
#' @slot metadata Contains a \code{list} to store any kind of experiment-wide
#' data and the metaoptions.

#' @details The \code{ProteomicsExperiment} object is just a ProteinExperiment
#' object and a PeptideExperiment object together.
#'
#' The rows of the \code{ProteinExperiment} object represents proteins. The rows
#' of the \code{PeptideExperiment} object represents peptides.
#'
#' The columns of the \code{ProteomicsExperiment} object represent samples.
#' Samples are shared at both protein and peptide levels.
#'
#' Experiment-wide information can be stored in the \code{metadata} slot, which
#' is accessed with the \code{metadata} function. This contains a \code{list}
#' object in which each item is left to the discretion of the user. Some
#' possible examples could be: data of the experiment, author, machine used,
#' etc.
#'
#' ProteomicsExperiment options are stored in the \code{metadata} slot.
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
#' See \link{ProteomicsExperiment-constructor} for details.
#'
#' @section Accessors:
#' See \link{ProteomicsExperiment-accessors} for details.
#'
#' @seealso \code{\link{ProteomicsExperiment-constructor}},
#'          \code{\link{ProteomicsExperiment-accessors}},
#'          \code{\link{ProteinExperiment}},
#'          \code{\link{PeptideExperiment}}

#' @export
.ProteomicsExperiment <- setClass('ProteomicsExperiment',

  slots = representation(ProteinExperiment = 'ProteinExperiment',
                         PeptideExperiment = 'PeptideExperiment',
                         colData = 'DataFrame',
                         linkerDf = 'data.frame',
                         metadata = 'list')

)


.valid.ProteomicsExperiment.ncol <- function(x) {

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

.valid.ProteomicsExperiment.colData <- function(x) {

  if (any(duplicated(colData(x)))) {
    return('colData cannot have duplicated entries')
  } else {
    return(NULL)
  }

}

## Wrapper for all the validity check functions
.valid.ProteomicsExperiment <- function(x) {

  c(.valid.ProteomicsExperiment.ncol(x),
    .valid.ProteomicsExperiment.metaoptions(x),
    .valid.ProteomicsExperiment.linkerDf(x),
    .valid.ProteomicsExperiment.colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteomicsExperiment', .valid.ProteomicsExperiment)

