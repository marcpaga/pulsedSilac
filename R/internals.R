###### Logical checks ==========================================================

##### hasAssays
setGeneric('hasRowData', function(x){
  standardGeneric('hasRowData')
})

setMethod('hasRowData', 'ProteinExperiment', function(x){

  rd <- rowData(x)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

})

setMethod('hasRowData', 'PeptideExperiment', function(x){

  rd <- rowData(x)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }

})

setMethod('hasRowData', 'ProteomicsExperiment', function(x){

  outVec <- c(NA, NA)
  rd <- rowData(x@ProteinExperiment)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    outVec[1] <- TRUE
  } else {
    outVec[1] <- FALSE
  }

  rd <- rowData(x@PeptideExperiment)
  if (is(rd, 'DataFrame') & ncol(rd) > 0) {
    outVec[2] <- TRUE
  } else {
    outVec[2] <- FALSE
  }

  return(outVec)
})


###### Metaoptions =============================================================
setGeneric('metaoptions', function(x){
  standardGeneric('metaoptions')
})

setMethod('metaoptions', 'ProteinExperiment', function(x){

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  return(metadata(x)[metaoptions_names])

})

setMethod('metaoptions', 'PeptideExperiment', function(x){

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  return(metadata(x)[metaoptions_names])

})

setMethod('metaoptions', 'ProteomicsExperiment', function(x){

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol',
                         'idColProt',
                         'idColPep',
                         'linkedSubset',
                         'subsetMode')

  return(metadata(x)[metaoptions_names])

})

###### Specific metaoptions ====================================================

hasMetaoption <- function(x, option) {

  op <- metadata(x)[[option]]
  if (is.na(op)) {
    txt <- sprintf('Not defined in metadata: %s.', option)
    stop(txt)
  } else {
    return(TRUE)
  }

}

metaoptionInColData <- function(x, option) {

  if (hasMetaoption(x, option)) {
    op <- metadata(x)[[option]]

    if (op %in% names(colData(x))) {
      return(TRUE)
    } else {
      txt <- sprintf('Column not found in colData: %s.', op)
      stop(txt)
    }
  }

}

giveMetaoption <- function(x, option) {

  if (metaoptionInColData(x, option)) {
    op <- metadata(x)[[option]]
    return(op)
  }

}
