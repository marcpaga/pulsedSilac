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

  if (any(c('conditionCol', 'timeCol',
            'replicateIntCol', 'replicateTimeCol') %in% names(metadata))) {
    stop('Metadata elements cannot have the following names: "conditionCol", ',
         '"timeCol", "replicateIntCol" or "replicateTimeCol"')
  }

  if (is.null(metadata)) {
    metadata <- metaoptions
  } else {
    metadata <- c(metadata, metaoptions)
  }



  se <- SummarizedExperiment(assays = assays,
                             colData = colData,
                             rowData = rowData,
                             metadata = metadata)

  return(.ProteinExperiment(se))

}
