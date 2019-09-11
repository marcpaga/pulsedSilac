#' @title SilacProteomicsExperiment constructor
#' @name SilacProteomicsExperiment-constructor
#'
#' @description Constructor function for the ProteomicsExperiment class object.
#' It requires at minimum a \code{SilacProteinExperiment} and a
#' \code{SilacPeptideExperiment}. If the colData, metadata and metaoptions have
#' been already defined in those it is not necessary to give them again.
#'
#' @param SilacProteinExperiment A \code{SilacProteinExperiment} object.
#' @param SilacPeptideExperiment A \code{SilacPeptideExperiment} object.
#' @param colData A \code{data.frame} with sample information like
#' conditions, replicates, etc. If not provided uses the colData slot from
#' the \code{SilacProteinExperiment} and \code{SilacPeptideExperiment}.
#' @param linkerDf A \code{data.frame} output from \code{\link{buildLinkerDf}}.
#' @param idColProt A \code{character} indicating which column from the
#' rowData (protein) should be used as ids. Should be the same used in
#' \code{\link{buildLinkerDf}}.
#' @param idColPept A \code{character} indicating which column from the
#' rowData (peptide) should be used as ids. Should be the same used in
#' \code{\link{buildLinkerDf}}.
#' @param linkedSubset A \code{logical} if subsetting should be linked between
#' proteins and peptide.
#' @param subsetMode A \code{character}, either 'protein' or 'peptide'
#' indicating which level should be used first when subsetting.
#' @param conditionCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment conditions.
#' @param timeCol A \code{character}, which indicates the column name
#' in colData(x) that defines the different experiment timepoints.
#' @param proteinCol A \code{character}, which indicates the column name
#' in rowData(x) that defines to which protein a peptide is assigned.
#' @param metadata A \code{list} to store any kind of experiment-wide
#' data; like authors, dates, machines used... If not provided uses the metadata
#' from the \code{SilacProteinExperiment} and \code{SilacPeptideExperiment}.
#'
#' @return An object of class \code{SilacProteomicsExperiment}.
#'
#' @section Class description:
#' See \link{SilacProteomicsExperiment-class} for details.
#'
#' @section Accessors:
#' See \link{SilacProteomicsExperiment-accessors} for details.
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
#' protExp <- SilacProteinExperiment(assays = assays_protein,
#'                                   rowData = rowData_protein,
#'                                   colData = colData,
#'                                   conditionCol = 'condition',
#'                                   timeCol = 'time')
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
#' ## list with the relationships
#' protein_to_peptide <- list(A = c('a', 'b'), B = c('c'), C = c('d', 'e'))
#' ## function to build the data.frame
#' linkerDf <- buildLinkerDf(protIDs = LETTERS[1:3],
#'                           pepIDs  = letters[1:5],
#'                           protToPep = protein_to_peptide)
#'
#' ProteomicsExp <- SilacProteomicsExperiment(SilacProteinExperiment = protExp,
#'                                            SilacPeptideExperiment = peptExp,
#'                                            linkerDf = linkerDf)
#'
#' @importFrom taRifx merge.list
#' @import methods S4Vectors
#' @export
SilacProteomicsExperiment <- function(SilacProteinExperiment,
                                 SilacPeptideExperiment,
                                 colData,
                                 linkerDf,
                                 metadata,
                                 idColProt = NA,
                                 idColPept = NA,
                                 linkedSubset = TRUE,
                                 subsetMode = 'protein',
                                 conditionCol = NA,
                                 timeCol = NA,
                                 proteinCol = NA) {

  ## initialize metaoptions list
  PEmetaoptions <- list(idColProt = idColProt,
                        idColPept = idColPept,
                        linkedSubset = linkedSubset,
                        subsetMode = subsetMode,
                        conditionCol = conditionCol,
                        timeCol = timeCol,
                        proteinCol = proteinCol)

  ## linkedSubset only relevant if linkerDf is there
  if (missing(linkerDf)) {
    PEmetaoptions[['linkedSubset']] <- FALSE
    linkerDf <- data.frame()
  }


  ## both experiments are mandatory, otherwise you would use the other classes
  if (missing(SilacProteinExperiment) | missing(SilacPeptideExperiment)) {

    stop('For a SilacProteomicsExperiment you need both a SilacProteinExperiment and a ',
         'SilacPeptideExperiment.')
  ## both present
  } else {

    ## if not given colData check that the two levels have the same colData and
    ## use the protein one
    if (missing(colData)) {
      colDataProt <- SilacProteinExperiment@colData
      colDataPept <- SilacPeptideExperiment@colData
      if (identical(colDataProt, colDataPept)) {
        colData <- colDataProt
      } else {
        stop('"colData" from SilacProteinExperiment and SilacPeptideExperiment are not',
             ' equal')
      }
    }

    ## if not given metadata check that the two levels have the same metadata
    ## and use the protein one
    if (missing(metadata)) {
      metadataProt <- SilacProteinExperiment@metadata
      metadataPept <- SilacPeptideExperiment@metadata

      metadata <- merge.list(metadataProt, metadataPept)
      metadata <- metadata[-which(names(metadata) %in% names(PEmetaoptions))]
    }
  }

  if (is(colData, 'data.frame')) {
    colData <- DataFrame(colData)
  }

  metaoptions_input <- mergeMetaoptions(PEmetaoptions,
                                        metaoptions(SilacProteinExperiment))
  metaoptions_input <- mergeMetaoptions(metaoptions_input,
                                        metaoptions(SilacProteinExperiment))

  metadata <- c(metadata, metaoptions_input)

  PE <- new(Class = 'SilacProteomicsExperiment',
            SilacProteinExperiment = SilacProteinExperiment,
            SilacPeptideExperiment = SilacPeptideExperiment,
            colData = colData,
            linkerDf = linkerDf,
            metadata = metadata)

  PE <- synchronizeMetaoptions(PE)

  return(PE)
}

