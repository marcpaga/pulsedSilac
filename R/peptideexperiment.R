#' @title PeptideExperiment constructor
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
                              metadata = NULL) {

  ## initialize the metadata
  metaoptions <- list(conditionCol = conditionCol,
                      timeCol = timeCol,
                      replicateIntCol = replicateIntCol,
                      replicateTimeCol = replicateTimeCol)

  if (is.null(metadata)) {
    metadata <- metaoptions
  } else {
    metadata <- c(metadata, metaoptions)
  }

  se <- SummarizedExperiment(assays = assays,
                             colData = colData,
                             rowData = rowData,
                             metadata = metadata)

  return(.PeptideExperiment(se))

}
