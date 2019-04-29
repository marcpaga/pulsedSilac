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




###### ProteinExperiment =======================================================



###### PeptideExperiment =======================================================



###### ProteomicsExperiment ====================================================


