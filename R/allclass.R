###### PROTEIN EXPERIMENT ======================================================

#' @rdname ProteinExperiment-class
#' @title  ProteinExperiment class
#'
#' @description S4 class that extends the \code{SummarizedExperiment} class.
#' This class is designed for proteomics data, more especifically protein level
#' data.
#'
#' @slot assays Contains a \code{list} of matrices (assays) with protein level
#' data.
#' @slot elementMetadata Contains a \code{data.frame} with protein feature data
#' like protein names, molecular weight, etc.
#' @slot colData Contains a \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @slot metadata Contains a \code{list} to store any kind of experiment-wide
#' data.
#' @slot metaoptions Contains a \code{list} to store configuration variables for
#' data processing, analysis and plotting purposes.
#' @slot NAMES Unused. Inherited from the SummarizedExperiment class.

#' @details The \code{ProteinExperiment} class is extremely similar to the
#' \code{SummarizedExperiment} class. The only addition is the metaoptions slot.
#' This class has been designed for use in this package and also in combination
#' with the other two classes from the package: the
#' \code{\link{PeptideExperiment}} and \code{\link{ProteomicsExperiment}}
#' classes.
#'
#' The rows of the \code{ProteinExperiment} object represent proteins. These
#' can be accessed with the \code{\link{rowData}} function.
#'
#' The columns of the \code{ProteinExperiment} object represent samples.
#' Informationabout the samples is stored in a \code{\link{DataFrame}} that
#' can be accessed with the \code{colData} function.
#'
#' Quantification data is stored in the \code{assays} slot, which us accessed
#' with the \code{assays} function. It contains a \code{list} of matrices that
#' must have the same number of columns and rows as the \code{ProteinExperiment}
#' object.
#'
#' Experiment-wide information can be stored in the \code{metadata} slot, which
#' is accessed with the \code{metadata} function. This contains a \code{list}
#' object in which each item is left to the discretion of the user. Some
#' possible examples could be: data of the experiment, author, machine used,
#' etc.
#'
#' ProteinExperiment options are stored in the \code{metaoptions} slot, which
#' is accessed through the \code{metatopions} function. This contains a
#' \code{list} with some parameters that are automatically initialized by the
#' constructor. Some parameters are mandatory for certain functions or
#' operations. The user can add or remove items at their discretion. These
#' parameters are meant to help automate certain pipeline or data analysis
#' steps. These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'   \item{replicateIntCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered expression replicates.}
#'   \item{replicateTimeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered time replicates.}
#'}
#'
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
#' @title  PeptideExperiment class
#'
#' @description S4 class that extends the \code{SummarizedExperiment} class.
#' This class is designed for proteomics data, more especifically peptide level
#' data.
#'
#' @slot assays Contains a \code{list} of matrices (assays) with peptide level
#' data.
#' @slot elementMetadata Contains a \code{data.frame} with peptide feature data
#' like protein names, molecular weight, etc.
#' @slot colData Contains a \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @slot metadata Contains a \code{list} to store any kind of experiment-wide
#' data.
#' @slot metaoptions Contains a \code{list} to store configuration variables for
#' data processing, analysis and plotting purposes.
#' @slot NAMES Unused. Inherited from the SummarizedExperiment class.

#' @details The \code{PeptideExperiment} class is extremely similar to the
#' \code{SummarizedExperiment} class. The only addition is the metaoptions slot.
#' This class has been designed for use in this package and also in combination
#' with the other two classes from the package: the
#' \code{\link{ProteinExperiment}} and \code{\link{ProteomicsExperiment}}
#' classes.
#'
#' The rows of the \code{PeptideExperiment} object represent peptides. These
#' can be accessed with the \code{\link{rowData}} function.
#'
#' The columns of the \code{PeptideExperiment} object represent samples.
#' Informationabout the samples is stored in a \code{\link{DataFrame}} that
#' can be accessed with the \code{colData} function.
#'
#' Quantification data is stored in the \code{assays} slot, which us accessed
#' with the \code{assays} function. It contains a \code{list} of matrices that
#' must have the same number of columns and rows as the \code{PeptideExperiment}
#' object.
#'
#' Experiment-wide information can be stored in the \code{metadata} slot, which
#' is accessed with the \code{metadata} function. This contains a \code{list}
#' object in which each item is left to the discretion of the user. Some
#' possible examples could be: data of the experiment, author, machine used,
#' etc.
#'
#' PeptideExperiment options are stored in the \code{metaoptions} slot, which
#' is accessed through the \code{metatopions} function. This contains a
#' \code{list} with some parameters that are automatically initialized by the
#' constructor. Some parameters are mandatory for certain functions or
#' operations. The user can add or remove items at their discretion. These
#' parameters are meant to help automate certain pipeline or data analysis
#' steps. These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'   \item{replicateIntCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered expression replicates.}
#'   \item{replicateTimeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered time replicates.}
#'   \item{proteinCol}{\code{character} indicating the column name of
#'   \code{rowData} that defines to which protein a peptide is assigned.}
#'}
#'
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @import SummarizedExperiment
#' @export
.PeptideExperiment <- setClass(Class = 'PeptideExperiment',
                               contains = 'ProteinExperiment'
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
#' @title  ProteomicsExperiment class
#'
#' @description S4 class that contains a ProteinExperiment object and
#' a PeptideExperiment object. The two objects are linked by a
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
#' data.
#' @slot metaoptions Contains a \code{list} to store configuration variables for
#' data processing, analysis and plotting purposes.

#' @details The \code{ProteomicsExperiment} class is just a ProteinExperiment
#' object and a PeptideExperiment object together.
#'
#' The rows of the \code{ProteinExperiment} object represents proteins. The rows
#' of the \code{PeptideExperiment} object represents peptides.
#'
#' The columns of the \code{ProteomicsExperiment} object represent samples.
#' Samples are shared at both protein and peptide levels.
#'
#'
#' Experiment-wide information can be stored in the \code{metadata} slot, which
#' is accessed with the \code{metadata} function. This contains a \code{list}
#' object in which each item is left to the discretion of the user. Some
#' possible examples could be: data of the experiment, author, machine used,
#' etc.
#'
#' ProteomicsExperiment options are stored in the \code{metaoptions} slot, which
#' is accessed through the \code{metatopions} function. This contains a
#' \code{list} with some parameters that are automatically initialized by the
#' constructor. Some parameters are mandatory for certain functions or
#' operations. The user can add or remove items at their discretion. These
#' parameters are meant to help automate certain pipeline or data analysis
#' steps. Here there are the  same metaoptions as in \code{ProteinExperiment}
#' and \code{PeptideExperiment} objects plus some extras.
#' These metaoptions are:
#' \describe{
#'   \item{conditionCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different experiment conditions.}
#'   \item{timeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines the different timepoints of the experiment.}
#'   \item{replicateIntCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered expression replicates.}
#'   \item{replicateTimeCol}{\code{character} indicating the column name of
#'   \code{colData(x)} that defines which samples of a condition can be
#'   considered time replicates.}
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
#' @export
.ProteomicsExperiment <- setClass('ProteomicsExperiment',

  slots = representation(ProteinExperiment = 'ProteinExperiment',
                         PeptideExperiment = 'PeptideExperiment',
                         colData = 'DataFrame',
                         linkerDf = 'data.frame',
                         metadata = 'list',
                         metaoptions = 'list')

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

