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
                                 replicateTimeCol = NA) {

  ## initialize metaoptions list
  metaoptions <- list(idColProt = idColProt,
                      idColPep = idColPep,
                      linkedSubset = linkedSubset,
                      subsetMode = subsetMode,
                      conditionCol = conditionCol,
                      timeCol = timeCol,
                      replicateIntCol = replicateIntCol,
                      replicateTimeCol = replicateTimeCol)

  ## linkedSubset only relevant if both levels are present
  if (any(missing(ProteinExperiment),
          missing(PeptideExperiment),
          missing(linkerDf))) {
    metaoptions[['linkedSubset']] <- FALSE
  }


  ## how to deal with the slots if only one level is present
  ## no protein data
  if (missing(ProteinExperiment)) {

    metaoptions[['subsetMode']] <- 'peptide'

    ## use the peptide data if not given
    if (missing(colData)) {
      colData <- PeptideExperiment@colData
    }
    if (missing(metadata)) {
      metadata <- PeptideExperiment@metadata
    }

  ## no peptide data
  } else if (missing(PeptideExperiment)){

    metaoptions[['subsetMode']] <- 'protein'

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
      if (identical(metadataProt, metadataPept)) {
        metadata <- metadataProt
      } else {
        stop('"metadata" from ProteinExperiment and PeptideExperiment are not',
             ' equal')
      }
    }
  }

  if (is.null(metadata)) {
    metadata <- metaoptions
  } else {
    metadata <- merge.list(metadata, metaoptions)
  }

  if (missing(linkerDf)) {
    linkerDf <- data.frame()
  }

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = ProteinExperiment,
            PeptideExperiment = PeptideExperiment,
            colData = colData,
            linkerDf = linkerDf,
            metadata = metadata)

  return(PE)
}

