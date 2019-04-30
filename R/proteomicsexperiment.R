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
                                 idColPep = NA,
                                 linkedSubset = TRUE,
                                 subsetMode = 'protein',
                                 conditionCol = NA,
                                 timeCol = NA,
                                 replicateIntCol = NA,
                                 replicateTimeCol = NA,
                                 proteinCol = NA) {

  ## initialize metaoptions list
  PEmetaoptions <- list(idColProt = idColProt,
                        idColPep = idColPep,
                        linkedSubset = linkedSubset,
                        subsetMode = subsetMode,
                        conditionCol = conditionCol,
                        timeCol = timeCol,
                        replicateIntCol = replicateIntCol,
                        replicateTimeCol = replicateTimeCol,
                        proteinCol = proteinCol)

  ## linkedSubset only relevant if both levels are present
  if (any(missing(ProteinExperiment),
          missing(PeptideExperiment),
          missing(linkerDf))) {
    PEmetaoptions[['linkedSubset']] <- FALSE
  }


  ## how to deal with the slots if only one level is present
  ## no protein data
  if (missing(ProteinExperiment)) {

    PEmetaoptions[['subsetMode']] <- 'peptide'

    ## use the peptide data if not given
    if (missing(colData)) {
      colData <- PeptideExperiment@colData
    }
    if (missing(metadata)) {
      metadata <- PeptideExperiment@metadata
    }

  ## no peptide data
  } else if (missing(PeptideExperiment)){

    PEmetaoptions[['subsetMode']] <- 'protein'

    ## use the protein data if not given
    if (missing(colData)) {
      colData <- ProteinExperiment@colData
    }
    if (missing(metadata)) {
      metadata <- ProteinExperiment@metadata
    }

  ## both levels of data given
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

  if (missing(linkerDf)) {
    linkerDf <- data.frame()
  }

  metaoptions <- .mergeMetaoptions(PEmetaoptions, ProteinExperiment@metaoptions)
  metaoptions <- .mergeMetaoptions(metaoptions, PeptideExperiment@metaoptions)

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = ProteinExperiment,
            PeptideExperiment = PeptideExperiment,
            colData = colData,
            linkerDf = linkerDf,
            metadata = metadata,
            metaoptions = metaoptions)

  return(PE)
}

