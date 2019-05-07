#' @title ProteomicsExperiment constructor
#'
#' @export
#' @importFrom taRifx merge.list
#' @import methods S4Vectors
ProteomicsExperiment <- function(ProteinExperiment,
                                 PeptideExperiment,
                                 colData,
                                 linkerDf,
                                 metadata,
                                 idColProt = NA,
                                 idColPept = NA,
                                 linkedSubset = TRUE,
                                 subsetMode = 'protein',
                                 conditionCol = NA,
                                 timeCol = NA,
                                 replicateIntCol = NA,
                                 replicateTimeCol = NA,
                                 proteinCol = NA) {

  ## initialize metaoptions list
  PEmetaoptions <- list(idColProt = idColProt,
                        idColPept = idColPept,
                        linkedSubset = linkedSubset,
                        subsetMode = subsetMode,
                        conditionCol = conditionCol,
                        timeCol = timeCol,
                        replicateIntCol = replicateIntCol,
                        replicateTimeCol = replicateTimeCol,
                        proteinCol = proteinCol)

  ## linkedSubset only relevant if linkerDf is there
  if (missing(linkerDf)) {
    PEmetaoptions[['linkedSubset']] <- FALSE
    linkerDf <- data.frame()
  }


  ## both experiments are mandatory, otherwise you would use the other classes
  if (missing(ProteinExperiment) | missing(PeptideExperiment)) {

    stop('For a ProteomicsExperiment you need both a ProteinExperiment and a ',
         'PeptideExperiment.')
  ## both present
  } else {

    ## if not given colData check that the two levels have the same colData and
    ## use the protein one
    if (missing(colData)) {
      colDataProt <- ProteinExperiment@colData
      colDataPept <- PeptideExperiment@colData
      if (identical(colDataProt, colDataPept)) {
        colData <- colDataProt
      } else {
        stop('"colData" from ProteinExperiment and PeptideExperiment are not',
             ' equal')
      }
    }

    ## if not given metadata check that the two levels have the same metadata
    ## and use the protein one
    if (missing(metadata)) {
      metadataProt <- ProteinExperiment@metadata
      metadataPept <- PeptideExperiment@metadata

      metadata <- merge.list(metadataProt, metadataPept)
    }
  }

  if (is(colData, 'data.frame')) {
    colData <- DataFrame(colData)
  }

  metaoptions <- mergeMetaoptions(PEmetaoptions, ProteinExperiment@metaoptions)
  metaoptions <- mergeMetaoptions(metaoptions, PeptideExperiment@metaoptions)

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = ProteinExperiment,
            PeptideExperiment = PeptideExperiment,
            colData = colData,
            linkerDf = linkerDf,
            metadata = metadata,
            metaoptions = metaoptions)

  PE <- synchronizeMetaoptions(PE)

  return(PE)
}

