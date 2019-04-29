#' Show method for ProteomicsExperiment
#'
#' @param object ProteomicsExperiment object
#' @importMethodsFrom SummarizedExperiment show
#' @keywords internal
setMethod('show',
          'ProteomicsExperiment',
          function(object){

  # someVals <- function(chars, max = 6) {
  #
  #   if (length(chars) > max) {
  #
  #     ## even number
  #     if ((max %% 2) == 0) {
  #
  #       outchars <- rep(NA, max + 1)
  #       ## put 3 dots in the middle
  #       outchars[max%/%2 + 1] <- '...'
  #
  #       ## fill the last 2 and first 2 positions
  #       outchars[c(1,2)] <- chars[c(1,2)]
  #       outchars[length(outchars)] <- chars[length(chars)]
  #       outchars[length(outchars)-1] <- chars[length(chars)-1]
  #
  #       ## if there are still NAs
  #       if (sum(is.na(outchars)) > 0) {
  #
  #         ## add values based on their fraction position
  #         ## calculate their fraction position
  #         posdec <- which(is.na(outchars))/length(outchars)
  #         charsdec <- seq_len(length(chars))/length(chars)
  #
  #         ## get which values are closest to the fraction postions
  #         charspos <- vapply(posdec, FUN.VALUE =  c(1L), FUN = function(x){
  #           return(which(abs(x-charsdec) == min(abs(x-charsdec)))[1])
  #         })
  #
  #         ## put them in the rigth position
  #         nas <- sum(is.na(outchars))
  #         naspos <- which(is.na(outchars))
  #         for (i in seq_len(nas)) {
  #           outchars[naspos[i]] <- chars[charspos][i]
  #         }
  #
  #         return(outchars)
  #       }
  #
  #       ## no NAs left
  #       return(outchars)
  #
  #     } else {
  #
  #       outchars <- rep(NA, max)
  #       ## put 3 dots in the middle
  #       outchars[max%/%2 + 1] <- '...'
  #
  #       ## fill the last 2 and first 2 positions
  #       outchars[c(1,2)] <- chars[c(1,2)]
  #       outchars[length(outchars)] <- chars[length(chars)]
  #       outchars[length(outchars)-1] <- chars[length(chars)-1]
  #
  #       ## if there are still NAs
  #       if (sum(is.na(outchars)) > 0) {
  #
  #         ## add values based on their fraction position
  #         ## calculate their fraction position
  #         posdec <- which(is.na(outchars))/length(outchars)
  #         charsdec <- seq_len(length(chars))/length(chars)
  #
  #         ## get which values are closest to the fraction postions
  #         charspos <- vapply(posdec, FUN.VALUE = c(1L), FUN = function(x){
  #           return(which(abs(x-charsdec) == min(abs(x-charsdec)))[1])
  #         })
  #
  #         ## put them in the rigth position
  #         nas <- sum(is.na(outchars))
  #         naspos <- which(is.na(outchars))
  #         for (i in seq_len(nas)) {
  #           outchars[naspos[i]] <- chars[charspos][i]
  #         }
  #
  #         return(outchars)
  #       }
  #
  #       ## no NAs left
  #       return(outchars)
  #
  #     }
  #   } else {
  #     return(chars)
  #   }
  #
  # }
  #
  # scat <- function(fmt, vals = character(), exdent = 2, ...) {
  #
  #   vals <- ifelse(nzchar(vals), vals, "''")
  #   lbls <- paste(someVals(vals), collapse=" ")
  #   txt <- sprintf(fmt, length(vals), lbls)
  #
  #   outtext <- strwrap(txt, exdent=exdent, ...)
  #   return(outtext)
  # }
  #
  # cat("Class:", class(object), "\n")
  #
  # # base numbers
  # cat('Feature numbers:','\n')
  #
  # if (hasAssays(object)[1]) {
  #   cat("  Proteins:", nrow(object@assays), '\n')
  # }
  #
  # if (hasAssays(object)[2]) {
  #   cat("  Peptides:", nrow(object@assaysPep), '\n')
  # }
  #
  # if (hasLinkerDf(object)) {
  #   cat("  Links:", nrow(object@linkerDf), '\n')
  # }
  #
  # if (hasRowDataProt(object)) {
  #   p <- scat('Protein features(%d): %s\n', colnames(rowDataProt(object)))
  #   cat(p, sep = '\n')
  # }
  #
  # if (hasRowDataPep(object)) {
  #   p <- scat('Peptide features(%d): %s\n', colnames(rowDataPep(object)))
  #   cat(p, sep = '\n')
  # }
  #
  # if (hasColData(object)) {
  #   if (ncol(colData(object)) > 0) {
  #
  #     ## experiment design
  #     cat('Experiment design:\n')
  #
  #     p <- scat("Samples: %s\n", nrow(colData(object)))
  #     cat(paste0('  ', p), sep = '\n')
  #
  #     for (i in seq_len(ncol(colData(object)))) {
  #       p <- scat(paste0(colnames(colData(object)[i]),"(%d): %s\n"),
  #                 as.character(colData(object)[,i]))
  #       cat(paste0('  ', p), sep = '\n')
  #     }
  #   }
  # }
  #
  # if (hasAssayNames(object)[1]) {
  #   cat(scat('Protein assays (%d): %s\n', assayNamesProt(object)), sep = '\n')
  # }
  #
  # if (hasAssayNames(object)[2]) {
  #   cat(scat('Peptide assays (%d): %s\n', assayNamesPep(object)), sep = '\n')
  # }
  #
  # if (length(metadata(object)) != 0) {
  #   cat(scat('Metadata (%d): %s\n', names(metadata(object))), sep = '\n')
  # }

  show(object@ProteinExperiment)
  show(object@PeptideExperiment)

})

