#' @rdname ProteinPeptideExperiment-accessors
#' @aliases subset,ProteinExperiment-method
#' @export
setMethod('subset', 'ProteinExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname ProteinPeptideExperiment-accessors
#' @aliases subset,PeptideExperiment-method
#' @export
setMethod('subset', 'PeptideExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases subset,ProteomicsExperiment-method
#' @export
setMethod('subset', 'ProteomicsExperiment', function(x, ...) {

  if (giveMetaoption(x, 'subsetMode') == 'protein') {
    return(subsetProt(x, ...))
  } else if (giveMetaoption(x, 'subsetMode') == 'peptide') {
    return(subsetPept(x, ...))
  }

})

###### subsetProt

#' @rdname ProteomicsExperiment-accessors
#' @name subsetProt
#' @aliases subsetProt,ProteinExperiment-method
#' @export
setMethod('subsetProt', 'ProteinExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname ProteomicsExperiment-accessors
#' @name subsetProt
#' @aliases subsetProt,ProteomicsExperiment-method
#' @export
setMethod('subsetProt', 'ProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@ProteinExperiment), ...))
  metaoptions(x)[['subsetMode']] <- 'protein'
  return(x[rows, ])

})

###### subsetPept

#' @rdname ProteomicsExperiment-accessors
#' @name subsetPept
#' @aliases subsetPept,PeptideExperiment-method
#' @export
setMethod('subsetPept', 'PeptideExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname ProteomicsExperiment-accessors
#' @name subsetPept
#' @aliases subsetPept,ProteomicsExperiment-method

#' @export
setMethod('subsetPept', 'ProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@PeptideExperiment), ...))
  metaoptions(x)[['subsetMode']] <- 'peptide'
  return(x[rows, ])

})

#' @rdname ProteomicsExperiment-accessors
#' @aliases [,ProteomicsExperiment-ANY-ANY-method
#' @export
setMethod('[', c('ProteomicsExperiment', 'ANY', 'ANY'),
          function(x, i, j, ..., drop = TRUE) {

  ## slots not affected by subset
  new.metadata <- metadata(x)
  new.metaoptions <- metaoptions(x)

  ## check the metaoptions -----------------------------------------------------
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
