#' Constructs a linkerDf that can be used as input when constructing a
#' SilacProteomicsExperiment
#'
#' Constructs a 4 column \code{data.frame} that contains the relationships
#' between proteins and peptides: which peptides belong to which proteins and
#' vice versa.
#'
#' This data frame is used in several functions and operations involving the
#' \code{SilacProteomicsExperiment} class. Especially in object merging and
#' subsetting. The arguments protIDs and pepIDs are mandatory, but only one of
#' the protToPep or pepToProt arguments is necessary to build the linkerDf.
#'
#' @param protIDs a \code{character} vector with unique ids that can be mapped
#' back to the proteins. Must be in the same order as in the rowDataProt data
#' frame.
#' @param pepIDs a \code{character} vector with unique ids that can be mapped
#' back to the peptides. Must be in the same order as in the rowDataPep data
#' frame.
#' @param protToPep a list with the same length as the number of protIDs. Every
#' entry of the list contains the peptide ids that are linked to that protein.
#' Items in the list must be in the same order as in protIDs.
#' @param pepToProt a list with the same length as the number of pepIDs Every
#' entry of the list contains the protein ids that are linked to that peptide.
#' Items in the list must be in the same order as in pepIDs
#'
#' @return A \code{data.frame} with the following 4 columns:
#' \describe{
#'   \item{protID}{Column with the protein IDs.}
#'   \item{pepID}{Column with the peptide IDs.}
#'   \item{protRow}{Column with row numbers of protein IDs.}
#'   \item{protID}{Column with the row numbers of peptide IDs.}
#' }
#'
#' @examples
#'
#' ## list with the relationships
#' protein_to_peptide <- list(A = c('a', 'b'), B = c('c'), C = c('d', 'e'))
#' ## function to build the data.frame
#' linkerDf <- buildLinkerDf(protIDs = LETTERS[1:3],
#'                          pepIDs  = letters[1:5],
#'                          protToPep = protein_to_peptide)
#' linkerDf
#'
#' @export
buildLinkerDf <- function(protIDs,
                          pepIDs,
                          protToPep,
                          pepToProt){


  if (missing(protIDs) | missing(pepIDs)) {
    stop('Both protIDs and pepIDs are necessary')
  }
  if (missing(protToPep) & missing(pepToProt)) {
    stop('You need at least one: protToPep or pepToProt')
  }

  if (!missing(protIDs) & !missing(pepIDs)) {

    if (!(length(protIDs) == length(unique(protIDs))) |
    !(length(pepIDs) == length(unique(pepIDs)))) {
      stop('Not all ids are unique!')
    }

    if (!missing(protToPep) & !missing(pepToProt)) {

      if (length(protIDs) != length(protToPep) |
          length(pepIDs) != length(pepToProt)) {
        stop('The ID vectors lengths do not match with the lists lengths')
      }

      # build two matrices and compare them

      linkerDfA <- matrix(data = c(rep(protIDs, vapply(protToPep, length, 1)),
                                   unname(unlist(protToPep))), ncol = 2)

      linkerDfA <- cbind(linkerDfA, match(linkerDfA[, 1], protIDs))
      linkerDfA <- cbind(linkerDfA, match(linkerDfA[, 2], pepIDs))

      colnames(linkerDfA) <- c('protID', 'pepID', 'protRow', 'pepRow')

      linkerDfA <- checkLinkerDf(linkerDfA)

      linkerDfA <- as.data.frame(linkerDfA)
      linkerDfA[, 3:4] <- apply(linkerDfA[, 3:4], 2, as.numeric)


      linkerDfB <- matrix(data = c(unname(unlist(pepToProt)),
                                   rep(pepIDs, vapply(pepToProt, length, 1))),
                          ncol = 2)

      linkerDfB <- cbind(linkerDfB, match(linkerDfB[, 1], protIDs))
      linkerDfB <- cbind(linkerDfB, match(linkerDfB[, 2], pepIDs))

      colnames(linkerDfB) <- c('protID', 'pepID', 'protRow', 'pepRow')

      linkerDfB <- checkLinkerDf(linkerDfB)

      linkerDfB <- as.data.frame(linkerDfB)
      linkerDfB[, 3:4] <- apply(linkerDfB[, 3:4], 2, as.numeric)

      # we use this method to compare them because they can just be disordered
      # but have the same relationships

      ch1 <- paste(linkerDfA[, 1], linkerDfA[, 2], sep = '-')
      ch2 <- paste(linkerDfB[, 1], linkerDfB[, 2], sep = '-')

      linkerDfA[, 1] <- as.character(linkerDfA[, 1])
      linkerDfA[, 2] <- as.character(linkerDfA[, 2])
      linkerDfA[, 3] <- as.numeric(linkerDfA[, 3])
      linkerDfA[, 4] <- as.numeric(linkerDfA[, 4])

      if (all(ch1 %in% ch2) & all(ch2 %in% ch1)) {
        return(linkerDfA)
      } else {
        stop('There are some inconsistencies with the data given')
      }

    } else if (!missing(protToPep)) {

      if (length(protIDs) != length(protToPep)) {
        stop('The ID vectors lengths do not match with the lists lengths')
      }

      # build the matrix
      linkerDf <- matrix(data = c(rep(protIDs, vapply(protToPep, length, 1)),
                                      unname(unlist(protToPep))),
                         ncol = 2)

      # add row data
      linkerDf <- cbind(linkerDf, match(linkerDf[,1], protIDs))
      linkerDf <- cbind(linkerDf, match(linkerDf[,2], pepIDs))

      colnames(linkerDf) <- c('protID', 'pepID', 'protRow', 'pepRow')

      linkerDf <- checkLinkerDf(linkerDf)

      linkerDf <- as.data.frame(linkerDf)
      linkerDf[,3:4] <- apply(linkerDf[,3:4], 2, as.numeric)

      linkerDf[, 1] <- as.character(linkerDf[, 1])
      linkerDf[, 2] <- as.character(linkerDf[, 2])
      linkerDf[, 3] <- as.numeric(linkerDf[, 3])
      linkerDf[, 4] <- as.numeric(linkerDf[, 4])

      return(linkerDf)

    } else if (!missing(pepToProt)) {

      if (length(pepIDs) != length(pepToProt)) {
        stop('The ID vectors lengths do not match with the lists lengths')
      }

      # build the matrix
      linkerDf <- matrix(data = c(unname(unlist(pepToProt)),
                                  rep(pepIDs, vapply(pepToProt, length, 1))),
                         ncol = 2)

      # add the row data
      linkerDf <- cbind(linkerDf, match(linkerDf[,1], protIDs))
      linkerDf <- cbind(linkerDf, match(linkerDf[,2], pepIDs))

      colnames(linkerDf) <- c('protID', 'pepID', 'protRow', 'pepRow')

      linkerDf <- checkLinkerDf(linkerDf)

      linkerDf <- as.data.frame(linkerDf)
      linkerDf[,3:4] <- apply(linkerDf[,3:4], 2, as.numeric)

      linkerDf[, 1] <- as.character(linkerDf[, 1])
      linkerDf[, 2] <- as.character(linkerDf[, 2])
      linkerDf[, 3] <- as.numeric(linkerDf[, 3])
      linkerDf[, 4] <- as.numeric(linkerDf[, 4])

      return(linkerDf)

    }

  }

}
