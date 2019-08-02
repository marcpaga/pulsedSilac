#' @title ProteinExperiment constructor
#'
#' @description Constructor function for the ProteinExperiment class object.
#' For more information about the class check
#' \link{ProteinExperiment-class}.
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
#' @param replicateIntCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different expression replicates.
#' @param replicateTimeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different tiempoint replicates.
#' @param metadata A \code{list} to store any kind of experiment-wide
#' data; like authors, dates, machines used...
#'
#' @return An object of class ProteinExperiment.
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
#' protExp <- ProteinExperiment(assays = assays_protein,
#'                              rowData = rowData_protein,
#'                              colData = colData,
#'                              conditionCol = 'condition',
#'                              timeCol = 'time')
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
ProteinExperiment <- function(assays,
                              rowData = NULL,
                              colData = NULL,
                              conditionCol = NA,
                              timeCol = NA,
                              replicateIntCol = NA,
                              replicateTimeCol = NA,
                              metadata = NULL) {

  ## initialize the metadata
  metaoptions <- list(conditionCol = conditionCol,
                      timeCol = timeCol,
                      replicateIntCol = replicateIntCol,
                      replicateTimeCol = replicateTimeCol)


  se <- SummarizedExperiment(assays = assays,
                             colData = colData,
                             rowData = rowData,
                             metadata = metadata)

  return(.ProteinExperiment(se, metaoptions = metaoptions))

}

# ProteinExperiment accessor for a ProteomicsExperiment object

#' @rdname classAccessors
#' @aliases ProtExp
#'
#' @usage NULL
#' @export
ProtExp <- function(x) {

  return(x@ProteinExperiment)

}
