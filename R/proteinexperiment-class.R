#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @export
.ProteinExperiment <- setClass(Class = 'ProteinExperiment',
                               contains = 'SummarizedExperiment'
)

###### VALIDITY ######
.valid.ProteinExperiment.metaoptions<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  ## missing one or more of the metaoptions
  if (!all(metaoptions_names %in% names(metadata(x)))) {
    return(paste('Incomplete metadata, one or more of the following is',
                 'missing: "conditionCol", "timeCol", "replicateIntCol",',
                 "replicateTimeCol"))
  }

  ## duplicated metaoptions entries
  if (sum(names(metadata(x)) %in% metaoptions_names) > 4) {
    return(paste('Duplicated entries for one or more of the folling:',
                 '"conditionCol", "timeCol", "replicateIntCol", "replicateTimeCol"'))
  }

  ## non-unique metaoptions
  metaoptions <- unlist(metadata(x)[metaoptions_names])
  if (all(is.na(metaoptions))) {
    return(NULL)
  }

  if (sum(is.na(metaoptions)) > 0) {
    metaoptions <- metaoptions[-which(is.na(metaoptions))]
  }
  if (any(duplicated(metaoptions))) {
    return('Two or more metaoptions target the same column of colData')
  }

  ## all validity is passed
  return(NULL)
}

.valid.ProteinExperiment.metaoptions_colData<- function(x) {

  metaoptions_names <- c('conditionCol',
                         'timeCol',
                         'replicateIntCol',
                         'replicateTimeCol')

  ## for each of the metaoptions check:
  ##   - if its NA then skip
  ##   - if its numeric check that it can find that column number in colData
  ##   - if its character check that colData has that column name
  ##   - if its another type then error

  for (metaoption in metaoptions_names) {
    mopt <- metadata(x)[[metaoption]]

    if (is.na(mopt)) {
      next
    }
    if (is.character(mopt)) {

      if (!mopt %in% colnames(colData(x))) {
        return(paste(metaoption, 'column does not match colData'))
      }

    } else if (is.numeric(mopt)) {

      if (ncol(colData(x)) < mopt) {
        return(paste(metaoption, 'column does not match colData'))
      }

    } else {
      return(paste(metaoption, 'must be character or numeric'))
    }

  }

  ## all validity is passed
  return(NULL)
}

## Wrapper for all the validity check functions
.valid.ProteinExperiment <- function(x) {

  c(.valid.ProteinExperiment.metaoptions(x),
    .valid.ProteinExperiment.metaoptions_colData(x))

}

#' @importFrom S4Vectors setValidity2
#' @keywords internal
setValidity2('ProteinExperiment', .valid.ProteinExperiment)
