#' @title SilacPeptideExperiment constructor
#' @name SilacPeptideExperiment-constructor
#'
#' @description Constructor function for the SilacPeptideExperiment class
#' object.
#'
#' @param assays A named \code{list} of matrices (assays) with peptide level
#' data.
#' @param rowData A \code{data.frame} with peptide feature data
#' like protein names, molecular weight, etc.
#' @param colData A \code{data.frame} with sample information like
#' conditions, replicates, etc.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param timeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment timepoints.
#' @param proteinCol A \code{character}, which indicates the column name
#' in rowData(x) that defines to which protein a peptide is assigned.
#' @param metadata A \code{list} to store any kind of experiment-wide
#' data; like authors, dates, machines used...
#'
#' @return An object of class \code{SilacPeptideExperiment}.
#'
#' @section Class description:
#' See \link{SilacPeptideExperiment-class} for details.
#'
#' @section Accessors:
#' See \link{SilacProteinPeptideExperiment-accessors} for details.
#'
#' @examples
#'
#' ## assays
#' assays_peptide <- list(expression = matrix(1:15, ncol = 3))
#'
#' ## colData
#' colData <- data.frame(sample = c('A1', 'A2', 'A3'),
#'                       condition = c('A', 'A', 'A'),
#'                       time = c(1, 2, 3))
#' ## rowData
#' rowData_peptide <- data.frame(pept_id = letters[1:5],
#'                               prot_id = c('A', 'A', 'B', 'C', 'C'))
#' ## construct the ProteinExperiment
#' peptExp <- SilacPeptideExperiment(assays = assays_peptide,
#'                                   rowData = rowData_peptide,
#'                                   colData = colData,
#'                                   conditionCol = 'condition',
#'                                   timeCol = 'time')
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
SilacPeptideExperiment <- function(assays,
                              rowData = NULL,
                              colData = NULL,
                              conditionCol = NA,
                              timeCol = NA,
                              proteinCol = NA,
                              metadata = NULL) {

  pe <- SilacProteinExperiment(assays = assays,
                          colData = colData,
                          rowData = rowData,
                          metadata = metadata,
                          conditionCol = conditionCol,
                          timeCol = timeCol)

  metadata(pe)[['proteinCol']] <- proteinCol
  pe <- .SilacPeptideExperiment(pe)

}

# SilacPeptideExperiment accessor for a SilacProteomicsExperiment object
# the @name ProteomicsExperiment-accessors is a bit hacky, but otherwise
# the doc of ProteomicsExperiment-accessors gets PeptExp as name...

#' @rdname SilacProteomicsExperiment-accessors
#' @name SilacProteomicsExperiment-accessors
#' @aliases PeptExp
#' @export
PeptExp <- function(x) {

  return(x@SilacPeptideExperiment)

}
