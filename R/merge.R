#' @export
setMethod('merge', 'ProteinExperiment', function(x, y,
                                                 by, by.x = by, by.y = by,
                                                 all = TRUE, ...) {

  ## argument checks -----
  ## by.x and by.y are mandatory, if missing try to guess them from rownames
  ## of rowData
  if (missing(by)) {
    by = intersect(colnames(rowData(x)), colnames(rowData(y)))
    if (is.null(by)) {
      stop('by undefined, which columns should be used for merging?')
    }
    if (missing(by.x)) {
      by.x <- by
    }
    if (missing(by.y)) {
      by.y <- by
    }
  }


  ### start merging ----

  ## colData should have the same columns, therefore just use rbind
  new.colData <- rbind(x@colData, y@colData)

  ## use dataframe method of merge, allows for some user customization
  rD.x <- as.data.frame(rowData(x))
  rD.y <- as.data.frame(rowData(y))
  new.rowData <- merge(x = rD.x, y = rD.y,
                       by.x = by.x, by.y = by.y, all = all, ...)

  ## if both experiments have names we merge them by name
  ## else we merge them in position order
  if (hasAssayNames(x) & hasAssayNames(y)) {
    ## merge assays based on names
    all_names <- unique(c(assayNames(x), assayNames(y)))
    for (i in seq_along(all_names)) {

      if (i == 1) {
        new.assays <- list()
        cols.x <- seq_len(ncol(x))
        cols.y <- seq(ncol(x) + 1, ncol(x) + ncol(y), 1)
      }
      ## based on the rowData merge apply it to the assays
      new.assayTemplate <- matrix(data = NA,
                                  ncol = ncol(x) + ncol(y),
                                  nrow = nrow(new.rowData))

      ## doing it separate for each experiment in case there are unique assays
      ## in each of them
      if (all_names[i] %in% assayNames(x)) {

        ## get the rows based on the old and new rowData
        rows <- match(rD.x[,by.x[1]], new.rowData[,by.x[1]])
        new.assayTemplate[rows, cols.x] <- assays(x)[[all_names[i]]]

      }

      if (all_names[i] %in% assayNames(y)) {

        rows <- match(rD.y[,by.y[1]], new.rowData[,by.y[1]])
        new.assayTemplate[rows, cols.y] <- assays(y)[[all_names[i]]]

      }

      new.assays[[i]] <- new.assayTemplate
      if (i == length(all_names)) {
        names(new.assays) <- all_names
      }


    }
  ## merge them by order
  } else {
    assay_num <- max(length(assays(x)), length(assays(y)))
    for (i in seq_len(assay_num)) {

      if (i == 1) {
        new.assays <- list()
        cols.x <- seq_len(ncol(x))
        cols.y <- seq(ncol(x) + 1, ncol(x) + ncol(y), 1)
      }
      ## based on the rowData merge apply it to the assays
      new.assayTemplate <- matrix(data = NA,
                                  ncol = ncol(x) + ncol(y),
                                  nrow = nrow(new.rowData))

      ## doing it separate for each experiment in case one experiment has more
      ## assays than the other one
      if (i <= length(assays(x))) {

        ## get the rows based on the old and new rowData
        rows <- match(rD.x[,by.x[1]], new.rowData[,by.x[1]])
        new.assayTemplate[rows, cols.x] <- assays(x)[[i]]

      }

      if (i <= length(assays(y))) {

        rows <- match(rD.y[,by.y[1]], new.rowData[,by.y[1]])
        new.assayTemplate[rows, cols.y] <- assays(y)[[i]]

      }

      new.assays[[i]] <- new.assayTemplate

    }
  }


  ## metaoptions and metadata is all from x
  metaopts <- metaoptions(x)
  PE <- ProteinExperiment(assays = new.assays,
                          rowData = new.rowData,
                          colData = new.colData,
                          conditionCol = metaopts[['conditionCol']],
                          timeCol = metaopts[['timeCol']],
                          replicateIntCol = metaopts[['replicateIntCol']],
                          replicateTimeCol = metaopts[['replicateTimeCol']],
                          metadata = x@metadata)

  return(PE)
})


#' @export
setMethod('merge', 'PeptideExperiment', function(x, y,
                                                 by, by.x = by, by.y = by,
                                                 all = TRUE, ...) {

  ## argument checks -----
  ## by.x and by.y are mandatory, if missing try to guess them from rownames
  ## of rowData
  if (missing(by)) {
    by = intersect(colnames(rowData(x)), colnames(rowData(y)))
    if (is.null(by)) {
      stop('by undefined, which columns should be used for merging?')
    }
    if (missing(by.x)) {
      by.x <- by
    }
    if (missing(by.y)) {
      by.y <- by
    }
  }


  ### start merging ----

  ## colData should have the same columns, therefore just use rbind
  new.colData <- rbind(x@colData, y@colData)

  ## use dataframe method of merge, allows for some user customization
  rD.x <- as.data.frame(rowData(x))
  rD.y <- as.data.frame(rowData(y))
  new.rowData <- merge(x = rD.x, y = rD.y,
                       by.x = by.x, by.y = by.y, all = all, ...)

  ## if both experiments have names we merge them by name
  ## else we merge them in position order
  if (hasAssayNames(x) & hasAssayNames(y)) {
    ## merge assays based on names
    all_names <- unique(c(assayNames(x), assayNames(y)))
    for (i in seq_along(all_names)) {

      if (i == 1) {
        new.assays <- list()
        cols.x <- seq_len(ncol(x))
        cols.y <- seq(ncol(x) + 1, ncol(x) + ncol(y), 1)
      }
      ## based on the rowData merge apply it to the assays
      new.assayTemplate <- matrix(data = NA,
                                  ncol = ncol(x) + ncol(y),
                                  nrow = nrow(new.rowData))

      ## doing it separate for each experiment in case there are unique assays
      ## in each of them
      if (all_names[i] %in% assayNames(x)) {

        ## get the rows based on the old and new rowData
        rows <- match(rD.x[,by.x[1]], new.rowData[,by.x[1]])
        new.assayTemplate[rows, cols.x] <- assays(x)[[all_names[i]]]

      }

      if (all_names[i] %in% assayNames(y)) {

        rows <- match(rD.y[,by.y[1]], new.rowData[,by.y[1]])
        new.assayTemplate[rows, cols.y] <- assays(y)[[all_names[i]]]

      }

      new.assays[[i]] <- new.assayTemplate
      if (i == length(all_names)) {
        names(new.assays) <- all_names
      }


    }
    ## merge them by order
  } else {
    assay_num <- max(length(assays(x)), length(assays(y)))
    for (i in seq_len(assay_num)) {

      if (i == 1) {
        new.assays <- list()
        cols.x <- seq_len(ncol(x))
        cols.y <- seq(ncol(x) + 1, ncol(x) + ncol(y), 1)
      }
      ## based on the rowData merge apply it to the assays
      new.assayTemplate <- matrix(data = NA,
                                  ncol = ncol(x) + ncol(y),
                                  nrow = nrow(new.rowData))

      ## doing it separate for each experiment in case one experiment has more
      ## assays than the other one
      if (i <= length(assays(x))) {

        ## get the rows based on the old and new rowData
        rows <- match(rD.x[,by.x[1]], new.rowData[,by.x[1]])
        new.assayTemplate[rows, cols.x] <- assays(x)[[i]]

      }

      if (i <= length(assays(y))) {

        rows <- match(rD.y[,by.y[1]], new.rowData[,by.y[1]])
        new.assayTemplate[rows, cols.y] <- assays(y)[[i]]

      }

      new.assays[[i]] <- new.assayTemplate

    }
  }


  ## metaoptions and metadata is all from x
  metaopts <- metaoptions(x)
  PE <- PeptideExperiment(assays = new.assays,
                          rowData = new.rowData,
                          colData = new.colData,
                          conditionCol = metaopts[['conditionCol']],
                          timeCol = metaopts[['timeCol']],
                          replicateIntCol = metaopts[['replicateIntCol']],
                          replicateTimeCol = metaopts[['replicateTimeCol']],
                          proteinCol = metaopts[['proteinCol']],
                          metadata = x@metadata)

  return(PE)


})

#' @export
setMethod('merge', 'ProteomicsExperiment', function(x, y,
                                                    by.prot,
                                                    by.prot.x = by.prot,
                                                    by.prot.y = by.prot,
                                                    by.pept,
                                                    by.pept.x = by.pept,
                                                    by.pept.y = by.pept,
                                                    all = TRUE, ...) {

  if (missing(by.prot)) {
    by.prot = intersect(colnames(rowDataProt(x)), colnames(rowDataProt(y)))
    if (is.null(by.prot)) {
      stop('by undefined, which columns should be used for merging?')
    }
    if (missing(by.prot.x)) {
      by.prot.x <- by.prot
    }
    if (missing(by.prot.y)) {
      by.prot.y <- by.prot
    }
  }

  if (missing(by.pept)) {
    by.pept = intersect(colnames(rowDataPept(x)), colnames(rowDataPept(y)))
    if (is.null(by.pept)) {
      stop('by undefined, which columns should be used for merging?')
    }
    if (missing(by.pept.x)) {
      by.pept.x <- by.pept
    }
    if (missing(by.pept.y)) {
      by.pept.y <- by.pept
    }
  }


  ## do each merge separately from their correspondent merge methods
  new.ProteinExperiment <- merge(x = x@ProteinExperiment,
                                 y = y@ProteinExperiment,
                                 by = by.prot,
                                 by.x = by.prot.x,
                                 by.y = by.prot.y,
                                 all = all, ...)

  new.PeptideExperiment <- merge(x = x@PeptideExperiment,
                                 y = y@PeptideExperiment,
                                 by = by.pept,
                                 by.x = by.pept.x,
                                 by.y = by.pept.y,
                                 all = all, ...)


  ## both have linkerDf
  if (nrow(x@linkerDf) != 0 & nrow(x@linkerDf) != 0) {
    new.linkerDf <- rbindLinkerDf(x = x@linkerDf,
                                  y = y@linkerDf)
  } else {
    new.linkerDf <- data.frame()
  }

  PE <- new(Class = 'ProteomicsExperiment',
            ProteinExperiment = new.ProteinExperiment,
            PeptideExperiment = new.PeptideExperiment,
            colData = colData(new.ProteinExperiment),
            linkerDf = new.linkerDf,
            metadata = x@metadata,
            metaoptions = x@metaoptions)

  return(PE)

})

