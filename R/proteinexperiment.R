#' @title ProteinExperiment constructor
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
