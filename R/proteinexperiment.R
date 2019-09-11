#' @title SilacProteinExperiment constructor
#' @name SilacProteinExperiment-constructor
#'
#' @description Constructor function for the SilacProteinExperiment class
#' object.
#'
#' @param assays A named \code{list} of matrices (assays) with protein level
#' data.
#' @param rowData A \code{data.frame} with protein feature data
#' like protein names, molecular weight, etc.
#' @param colData A \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param timeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment timepoints.
#' @param metadata A \code{list} to store any kind of experiment-wide
#' data; like authors, dates, machines used...
#'
#' @return An object of class \code{SilacProteinExperiment}.
#'
#' @section Class description:
#' See \link{SilacProteinExperiment-class} for details.
#'
#' @section Accessors:
#' See \link{SilacProteinPeptideExperiment-accessors} for details.
#'
#' @examples
#'
#' ## assays
#' assays_protein <- list(expression = matrix(1:9, ncol = 3))
#'
#' ## colData
#' colData <- data.frame(sample = c('A1', 'A2', 'A3'),
#'                       condition = c('A', 'A', 'A'),
#'                     time = c(1, 2, 3))
#' ## rowData
#' rowData_protein <- data.frame(prot_id = LETTERS[1:3])
#'
#' ## construct the ProteinExperiment
#' protExp <- SilacProteinExperiment(assays = assays_protein,
#'                                   rowData = rowData_protein,
#'                                   colData = colData,
#'                                   conditionCol = 'condition',
#'                                   timeCol = 'time')
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SilacProteinExperiment <- function(assays,
                              rowData = NULL,
                              colData = NULL,
                              conditionCol = NA,
                              timeCol = NA,
                              metadata = NULL) {

  metaoptions_names <- c('conditionCol', 'timeCol')
  if (any(metaoptions_names %in% names(metadata))) {
    stop('"conditionCol" and "timeCol" are metaoptions reserved names')
  }

  se <- SummarizedExperiment(assays = assays,
                             colData = colData,
                             rowData = rowData,
                             metadata = metadata)

  metadata(se)[['conditionCol']] <- conditionCol
  metadata(se)[['timeCol']] <- timeCol

  return(.SilacProteinExperiment(se))

}

# ProteinExperiment accessor for a ProteomicsExperiment object

#' @rdname SilacProteomicsExperiment-accessors
#'
#' @export
ProtExp <- function(x) {

  return(x@SilacProteinExperiment)

}
