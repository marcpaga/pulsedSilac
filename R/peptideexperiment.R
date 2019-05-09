#' @title PeptideExperiment constructor
#'
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
