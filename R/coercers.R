## Coercers to SummarizedExperiment objects ------------------------------------

# works for both ProteinExperiment and PeptideExperiment

#' @rdname ProteinPeptideExperiment-accessors
#' @name coerce
#' @aliases coerce,ProteinExperiment,SummarizedExperiment-method
#' @exportMethod coerce
setAs('ProteinExperiment', 'SummarizedExperiment',
      function(from){

        SummarizedExperiment(assays = assays(from),
                             rowData = rowData(from),
                             colData = colData(from),
                             metadata = metadata(from))

      }
)


## Coercers to data.frame objects ----------------------------------------------

# works for both ProteinExperiment and PeptideExperiment

#' @rdname ProteinPeptideExperiment-accessors
#' @name coerce
#' @aliases coerce,ProteinExperiment,data.frame-method
#' @exportMethod coerce
setAs('ProteinExperiment', 'data.frame',
      function(from){

        n_samples <- ncol(from)

        ## define sample names
        sampleNames <- rownames(colData(from))
        if (is.null(sampleNames)) {
          coldata <- vapply(colData(from), as.character,
                            rep('', n_samples))
          sampleNames <- apply(coldata, 1, function(x){
            paste(x, collapse = '.')
          })
        }

        ## define assay names
        namesAssay <- assayNames(from)
        if (is.null(namesAssay)) {
          namesAssay <- LETTERS[seq_len(length(assays(from)))]
        }

        ## merge the two to make column names for the final data.frame
        colNames <- vapply(namesAssay, FUN.VALUE = rep('', n_samples),
                           FUN = function(x){
                             paste(x, sampleNames, sep = '_')
                           })
        colNames <- as.vector(colNames)

        assaysMat <- do.call('cbind', as.list(assays(from)))
        colnames(assaysMat) <- colNames

        finalDf <- cbind(as.data.frame(rowData(from)),
                         as.data.frame(assaysMat))
        return(finalDf)
      }
)

#' @rdname ProteinPeptideExperiment-accessors
#' @name coerce
#' @aliases coerce,SummarizedExperiment,ProteinExperiment-method
#' @exportMethod coerce
setAs('ProteinExperiment', 'SummarizedExperiment',
      function(from){

        ProteinExperiment(assays = assays(from),
                          rowData = rowData(from),
                          colData = colData(from),
                          metadata = metadata(from))

      }
)
