#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases subset,SilacProteinExperiment-method
#' @export
setMethod('subset', 'SilacProteinExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname SilacProteinPeptideExperiment-accessors
#' @aliases subset,SilacPeptideExperiment-method
#' @export
setMethod('subset', 'SilacPeptideExperiment', function(x, ...) {

  rows <- which(with(rowData(x), ...))
  return(x[rows, ])

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases subset,SilacProteomicsExperiment-method
#' @export
setMethod('subset', 'SilacProteomicsExperiment', function(x, ...) {

  if (.giveMetaoption(x, 'subsetMode') == 'protein') {
    return(subsetProt(x, ...))
  } else if (.giveMetaoption(x, 'subsetMode') == 'peptide') {
    return(subsetPept(x, ...))
  }

})

###### subsetProt

#' @rdname SilacProteomicsExperiment-accessors
#' @name subsetProt
#' @aliases subsetProt,SilacProteinExperiment-method
#' @export
setMethod('subsetProt', 'SilacProteinExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name subsetProt
#' @aliases subsetProt,SilacProteomicsExperiment-method
#' @export
setMethod('subsetProt', 'SilacProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@SilacProteinExperiment), ...))
  metadata(x)[['subsetMode']] <- 'protein'
  return(x[rows, ])

})

###### subsetPept

#' @rdname SilacProteomicsExperiment-accessors
#' @name subsetPept
#' @aliases subsetPept,SilacPeptideExperiment-method
#' @export
setMethod('subsetPept', 'SilacPeptideExperiment', function(x, ...) {

  return(subset(x, ...))

})

#' @rdname SilacProteomicsExperiment-accessors
#' @name subsetPept
#' @aliases subsetPept,SilacProteomicsExperiment-method

#' @export
setMethod('subsetPept', 'SilacProteomicsExperiment', function(x, ...) {

  rows <- which(with(rowData(x@SilacPeptideExperiment), ...))
  metadata(x)[['subsetMode']] <- 'peptide'
  return(x[rows, ])

})

#' @rdname SilacProteomicsExperiment-accessors
#' @aliases [,SilacProteomicsExperiment-ANY-ANY-method
#' @export
setMethod('[', c('SilacProteomicsExperiment', 'ANY', 'ANY'),
          function(x, i, j, ..., drop = TRUE) {

  ## slots not affected by subset
  new.metadata <- metadata(x)

  ## check the metaoptions -----------------------------------------------------
  ## should linked subsetting be used
  linked <- .giveMetaoption(x, 'linkedSubset')

  if (linked) {
    if (nrow(linkerDf(x)) == 0) {
      warning('linkedDf not present, linkedSubset set to FALSE')
      linked <- FALSE
    }
  }

  if (linked) {
    idColProt <- try(.giveMetaoption(x, 'idColProt'))
    idColPept <- try(.giveMetaoption(x, 'idColPept'))

    if (is(idColProt, 'try-error') | is(idColProt, 'try-error')) {
      txt <- paste('idColProt and/or idColPep in metaoptions are not specified',
                   'the linkedSubset set to FALSE')
      warning(txt)
      linked <- FALSE
    }
  }

  ## subset based on protein or peptide
  subsetMode <- .giveMetaoption(x, 'subsetMode')

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
    new.ProteinExperiment <- x@SilacProteinExperiment[i, j]

    if (linked) {

      ## subset the linkerDf too
      new.linkerDf <- subsetLinkerDf(x = linkerDf(x), protRows = i)
      ## for the peptide part get which peptide rows belong to the proteins
      ## and use the SummarizedExperiment method too
      peptiderows <- unique(new.linkerDf[,'pepRow'])
      new.PeptideExperiment <- x@SilacPeptideExperiment[peptiderows, j]

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
      new.PeptideExperiment <- x@SilacPeptideExperiment[, j]
    }


    PE <- new(Class = 'SilacProteomicsExperiment',
              SilacProteinExperiment = new.ProteinExperiment,
              SilacPeptideExperiment = new.PeptideExperiment,
              colData = new.colData,
              linkerDf = new.linkerDf,
              metadata = new.metadata)

    return(PE)

  }


  if (subsetMode == 'peptide') {

    ## for the protein part use the SummarizedExperiment method
    new.PeptideExperiment <- x@SilacPeptideExperiment[i, j]

    if (linked) {

      ## subset the linkerDf too
      new.linkerDf <- subsetLinkerDf(x = linkerDf(x), pepRows = i)
      ## for the peptide part get which peptide rows belong to the proteins
      ## and use the SummarizedExperiment method too
      proteinrows <- unique(new.linkerDf[, 'protRow'])
      new.ProteinExperiment <- x@SilacProteinExperiment[proteinrows, j]

      ## recalculate the new rows for the linkerDf
      newProtRows <- match(new.linkerDf[,'protID'],
                           rowData(new.ProteinExperiment)[, idColProt])
      new.linkerDf[, 'protRow'] <- newProtRows

      pepids <- as.character(rowData(new.PeptideExperiment)[, idColPept])
      newPepRows <- match(new.linkerDf[, 'pepID'], pepids)
      new.linkerDf[, 'pepRow'] <- newPepRows

    } else {
      new.linkerDf <- linkerDf(x)
      new.ProteinExperiment <- x@SilacProteinExperiment[, j]
    }


    PE <- new(Class = 'SilacProteomicsExperiment',
              SilacProteinExperiment = new.ProteinExperiment,
              SilacPeptideExperiment = new.PeptideExperiment,
              colData = new.colData,
              linkerDf = new.linkerDf,
              metadata = new.metadata)

    return(PE)

  }

})
