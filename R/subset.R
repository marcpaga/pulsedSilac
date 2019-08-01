###### subset

#' @rdname subset
#' @name subset
#'
#' @title Subset functions
#'
#' @description Subset functions for the ProteinExperiment, PeptideExperiment
#' and ProteomicsExperiment classes.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param ... For \code{subset}, \code{subsetProt} and \code{subsetPep}
#' it is a logical comparison using a column name from the respective rowData
#' \code{data.frame}. Otherwise unused.
#' @param i,j For \code{`[`}, \code{i}, \code{j} are subscripts that can act to
#' subset the rows and columns of \code{x}.
#' @param drop Unused.
#'
#' @details For the ProteinExperiment and PeptideExperiment classes these work
#' in the same way as in the SummarizedExperiment class. For the subset
#' functions a logical comparison is done on the rowData \code{data.frame} and
#' the rows (proteins or peptides) that are TRUE are kept in the returned
#' object. For the `[` operator, numeric values indicating rows and columns
#' can be used to subset.
#'
#' The ProteomicsExperiment class is a bit more complex since there are two
#' levels at which the subset can be done and these two levels can be linked or
#' not.
#'
#' If the metaoption 'linkedSubset' is TRUE, then when subsetting on one level,
#' the proteins/peptide linked to such level are also subsetted. Otherwise, one
#' of the levels remains unmodified.
#'
#' \code{subsetProt} can be used to apply \code{subset} at the rowData
#' \code{data.frame} of the protein level. \code{subsetProt} can be used to
#' apply \code{subset} at the rowData \code{data.frame} of the peptide level. If
#' \code{subset} is used, then \code{subsetProt} or \code{subsetPept} will be
#' used depending on the metaoption 'subsetMode'.
#'
#' `[` acts in the same manner as calling \code{subset}. In this case numerics
#' are used and samples can also be selected.
#'
#' The vignette offers a detailed set of simple examples with all the
#' possible cases.
#'
#' @return A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#'
#' @examples
#' wormsPE[1,1]
#' subsetProt(wormsPE, protein_id == 'AC3.2')
#' subsetPept(wormsPE, Sequence == 'AIQEISDYHFLIK')
NULL

#' @rdname subset
#' @export
setMethod('subset', 'ProteinExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname subset
#' @export
setMethod('subset', 'PeptideExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname subset
#' @export
setMethod('subset', 'ProteomicsExperiment', function(x, ...) {

  if (giveMetaoption(x, 'subsetMode') == 'protein') {
    return(subsetProt(x, ...))
  } else if (giveMetaoption(x, 'subsetMode') == 'peptide') {
    return(subsetPept(x, ...))
  }

})

###### subsetProt

#' @rdname subset
#' @export
setMethod('subsetProt', 'ProteinExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname subset
#' @export
setMethod('subsetProt', 'ProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@ProteinExperiment), ...))
  metaoptions(x)[['subsetMode']] <- 'protein'
  return(x[rows, ])

})

###### subsetPept

#' @rdname subset
#' @export
setMethod('subsetPept', 'PeptideExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname subset
#' @export
setMethod('subsetPept', 'ProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@PeptideExperiment), ...))
  metaoptions(x)[['subsetMode']] <- 'peptide'
  return(x[rows, ])

})

#' @rdname subset
#' @export
setMethod('[', c('ProteinExperiment', 'ANY', 'ANY'),
          function(x, i, j, ..., drop = TRUE) {

  return(callNextMethod())

})

#' @rdname subset
#' @export
setMethod('[', c('PeptideExperiment', 'ANY', 'ANY'),
          function(x, i, j, ..., drop = TRUE) {

  return(callNextMethod())

})

#' @rdname subset
#' @export
setMethod('[', c('ProteomicsExperiment', 'ANY', 'ANY'),
          function(x, i, j, ..., drop = TRUE) {


  ## slots not affected by subset
  new.metadata <- metadata(x)
  new.metaoptions <- metaoptions(x)

  ## should linked subsetting be used
  linked <- giveMetaoption(x, 'linkedSubset')

  if (linked) {
    if (nrow(linkerDf(x)) == 0) {
      warning('linkedDf not present, linkedSubset set to FALSE')
      linked <- FALSE
    }
  }

  if (linked) {
    idColProt <- try(giveMetaoption(x, 'idColProt'))
    idColPept <- try(giveMetaoption(x, 'idColPept'))

    if (is(idColProt, 'try-error') | is(idColProt, 'try-error')) {
      txt <- paste('idColProt and/or idColPep in metaoptions are not specified',
                   'the linkedSubset set to FALSE')
      warning(txt)
      linked <- FALSE
    }
  }

  ## subset based on protein or peptide
  subsetMode <- giveMetaoption(x, 'subsetMode')

  ## in case of missing values
  if (missing(i)) {

    if (subsetMode == 'protein') {
      i <- seq_len(nrow(x)[1])

    } else if (subsetMode == 'peptide') {
      i <- seq_len(nrow(x)[2])
    }
  }

  if (missing(j)) {
    j <- seq_len(ncol(x))
  }

  ## colData is protein/peptide independent
  new.colData <- colData(x)[j, , drop = FALSE]


  if (subsetMode == 'protein') {

    ## for the protein part use the SummarizedExperiment method
    new.ProteinExperiment <- x@ProteinExperiment[i, j]

    if (linked) {

      ## subset the linkerDf too
      new.linkerDf <- subsetLinkerDf(x = linkerDf(x), protRows = i)
      ## for the peptide part get which peptide rows belong to the proteins
      ## and use the SummarizedExperiment method too
      peptiderows <- unique(new.linkerDf[,'pepRow'])
      new.PeptideExperiment <- x@PeptideExperiment[peptiderows, j]

      ## recalculate the new rows for the linkerDf
      ## protein part
      newProtRows <- match(new.linkerDf[,'protID'],
                           rowData(new.ProteinExperiment)[, idColProt])
      new.linkerDf[, 'protRow'] <- newProtRows

      ## peptide part
      pepids <- as.character(rowData(new.PeptideExperiment)[, idColPept])
      newPepRows <- match(new.linkerDf[, 'pepID'], pepids)
      new.linkerDf[, 'pepRow'] <- newPepRows

    } else {
      new.linkerDf <- linkerDf(x)
      new.PeptideExperiment <- x@PeptideExperiment[, j]
    }


    PE <- new(Class = 'ProteomicsExperiment',
              ProteinExperiment = new.ProteinExperiment,
              PeptideExperiment = new.PeptideExperiment,
              colData = new.colData,
              linkerDf = new.linkerDf,
              metadata = new.metadata,
              metaoptions = new.metaoptions)

    return(PE)

  }


  if (subsetMode == 'peptide') {

    ## for the protein part use the SummarizedExperiment method
    new.PeptideExperiment <- x@PeptideExperiment[i, j]

    if (linked) {

      ## subset the linkerDf too
      new.linkerDf <- subsetLinkerDf(x = linkerDf(x), pepRows = i)
      ## for the peptide part get which peptide rows belong to the proteins
      ## and use the SummarizedExperiment method too
      proteinrows <- unique(new.linkerDf[, 'protRow'])
      new.ProteinExperiment <- x@ProteinExperiment[proteinrows, j]

      ## recalculate the new rows for the linkerDf
      newProtRows <- match(new.linkerDf[,'protID'],
                           rowData(new.ProteinExperiment)[, idColProt])
      new.linkerDf[, 'protRow'] <- newProtRows

      pepids <- as.character(rowData(new.PeptideExperiment)[, idColPept])
      newPepRows <- match(new.linkerDf[, 'pepID'], pepids)
      new.linkerDf[, 'pepRow'] <- newPepRows

    } else {
      new.linkerDf <- linkerDf(x)
      new.ProteinExperiment <- x@ProteinExperiment[, j]
    }


    PE <- new(Class = 'ProteomicsExperiment',
              ProteinExperiment = new.ProteinExperiment,
              PeptideExperiment = new.PeptideExperiment,
              colData = new.colData,
              linkerDf = new.linkerDf,
              metadata = new.metadata,
              metaoptions = new.metaoptions)

    return(PE)

  }

})
