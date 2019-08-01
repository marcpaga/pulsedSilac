#' @title PeptideExperiment constructor
#'
#' @description Constructor function for the PeptideExperiment class object.
#' For more information about the class check
#' \link{PeptideExperiment-class}.
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
#' @param replicateIntCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different expression replicates.
#' @param replicateTimeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different tiempoint replicates.
#' @param proteinCol A \code{character}, which indicates the column name
#' in rowData(x) that defines to which protein a peptide is assigned.
#' @param metadata A \code{list} to store any kind of experiment-wide
#' data; like authors, dates, machines used...
#'
#' @return An object of class PeptideExperiment
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
#' peptExp <- PeptideExperiment(assays = assays_peptide,
#'                              rowData = rowData_peptide,
#'                              colData = colData,
#'                              conditionCol = 'condition',
#'                              timeCol = 'time')
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
PeptideExperiment <- function(assays,
                              rowData = NULL,
                              colData = NULL,
                              conditionCol = NA,
                              timeCol = NA,
                              replicateIntCol = NA,
                              replicateTimeCol = NA,
                              proteinCol = NA,
                              metadata = NULL) {

  ## initialize the metadata
  metaoptions <- list(conditionCol = conditionCol,
                      timeCol = timeCol,
                      replicateIntCol = replicateIntCol,
                      replicateTimeCol = replicateTimeCol,
                      proteinCol = proteinCol)


  pe <- ProteinExperiment(assays = assays,
                          colData = colData,
                          rowData = rowData,
                          metadata = metadata,
                          conditionCol = conditionCol,
                          timeCol = timeCol,
                          replicateIntCol = replicateIntCol,
                          replicateTimeCol = replicateTimeCol)

  metaoptions(pe)[['proteinCol']] <- proteinCol
  pe <- .PeptideExperiment(pe)

}
