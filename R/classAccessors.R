# Documentation for the getters and setters for the ProteinExperiment,
# PeptideExperiment and ProteomicsExperiment classes

#' @rdname classAccessors
#' @name classAccessors
#'
#' @title pulsedSilac classes accessors
#'
#' @description All the accessors, dimension, subsetting and merging functions
#' that work on \code{ProteinExperiment}, \code{PeptideExperiment} and
#' \code{ProteomicsExperiment} classes. For the first two, methods for
#' \code{SummarizedExperiment} should work too. Most functions work the same as
#' in \code{SummarizedExperiment} and are also shown in the vignette. For
#' \code{subset}, \code{metaoptions} and \code{merge} there are additional
#' documentation entries with a more detailed explanations.
#'
#' @param x A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param y A \code{ProteinExperiment}, \code{PeptideExperiment} or a
#' \code{ProteomicsExperiment} object.
#' @param i,j For \code{`[`}, \code{i}, \code{j} are subscripts that can act to
#' subset the rows and columns of \code{x}.
#' @param drop A \code{logical} indicating if dimensions should be lowered if
#' possible when subsetting.
#' @param ... For \code{rbind} and \code{cbind} are \code{ProteinExperiment},
#' \code{PeptideExperiment} or \code{ProteomicsExperiment} objects to be
#' joined together. For \code{subset}, \code{subsetProt} and \code{subsetPep}
#' it is a logical comparison using a column name from the respective rowData
#' \code{data.frame}. Otherwise unused.
#' @param value An object of class specified in the S4 method signature or as
#' described in the following sections.
#' @param name Column name of colData.
#' @param use.names Unused.
#' @param withDimnames Unused.
#' @param deparse.level Unused.
#'
#' @details For \code{ProteinExperiment} and \code{PeptideExperiment} most
#' functions use the same generic as in \code{SummarizedExperiment}. For
#' \code{ProteomicsExperiment} new generics with and added suffix (Prot or Pept)
#' to indicate which level should be used. These new generics also work
#' on the previous two classes.
#'
#' @usage
#' ## Accessors
#'
#' assays(x, ..., withDimnames = TRUE)
#' assays(x, ..., withDimnames = TRUE) <- value
#' assaysProt(x)
#' assaysProt(x) <- value
#' assaysPept(x)
#' assaysPept(x) <- value
#'
#' assayNames(x, ...)
#' assayNames(x, ...) <- value
#' assayNamesProt(x)
#' assayNamesProt(x) <- value
#' assayNamesPept(x)
#' assayNamesPept(x) <- value
#'
#' colData(x, ...)
#' colData(x, ...) <- value
#'
#' rowData(x, use.names = TRUE, ...)
#' rowData(x, ...) <- value
#' rowDataProt(x)
#' rowDataProt(x) <- value
#' rowDataPept(x)
#' rowDataPept(x) <- value
#'
#' metadata(x, ...)
#' metadata(x, ...) <- value
#'
#' #metaoptions(x)
#' #metaoptions(x) <- value
#'
#' linkerDf(x)
#' linkerDf(x) <- value
#'
#' ## Dimensions
#'
#' dim(x)
#' length(x)
#' nrow(x)
#' ncol(x)
#'
#' ## Subsetting
#'
#' x$name
#' x$name <- value
#'
#' #x[i, j, ..., drop = TRUE]
#' #subset(x, ...)
#' #subsetProt(x, ...)
#' #subsetPept(x, ...)
#'
#' ## Merging
#'
#' cbind(..., deparse.level = 1)
#' rbind(..., deparse.level = 1)
#' #merge(x, y, by, by.x = by, by.y = by, all = TRUE, ...)
#'
#'
#' @section Accessors:
#'
#' The following functions can be used to access the data in the class slots
#'
#' \describe{
#'   \item{\code{assays}, \code{assaysProt}, \code{assaysPept}:}{Access the
#'    assays (list of matrices) of the object. Value should be a matrix or
#'    list of matrices.}
#'   \item{\code{assayNames}, \code{assayNamesProt}, \code{assayNamesPept}:}{
#'   Access the assay names of the object. Value should be a character vector.}
#'   \item{\code{rowData}, \code{rowDataProt}, \code{rowDataPept}:}{Access the
#'   protein/peptide feature \code{data.frame} of the object. Value should be a
#'   \code{data.frame} with as many rows as proteins/peptides.}
#'   \item{\code{colData}:}{Access the samples \code{data.frame} of the object.
#'   Value should be a \code{data.frame} with as many rows as samples.}
#'   \item{\code{metadata}:}{Access the metadata \code{list} of the object.
#'   Value should be a \code{list}.}
#'   \item{\code{metaoptions}:}{Access the metaoptions \code{list} of the
#'   object. Value should be a \code{list}.}
#'   \item{\code{linkerDf}:}{Access the linker \code{data.frame} of the object
#'   (only ProteomicsExperiment). Value should be a \code{data.frame} output
#'   from \code{\link{buildLinkerDf}}.}
#' }
#'
#' @section Dimensions:
#'
#' The following functions can be used to get the number of proteins/peptides
#' and number of samples:
#'
#' \describe{
#'   \item{\code{nrow}:}{Gives how many proteins and/or peptides the object
#'   has.}
#'   \item{\code{ncol}:}{Gives how many sample the object has.}
#'   \item{\code{dim}:}{Gives both the number of proteins/peptides and the
#'   number of samples the object has. For \code{ProteinExperiment} and
#'   \code{PeptideExperiment} this is a numeric vector and for
#'   \code{ProteomicsExperiment} it is a numeric matrix.}
#'   \item{\code{length}:}{Gives how many proteins and/or peptides the object
#'   has.}
#' }
#'
#' @section Subsetting:
#'
#' The following functions can be used to subset the different classes:
#'
#' \describe{
#'   \item{\code{$}:}{Gives a column from colData by name.}
#'   \item{\code{`[`}:}{Can be used to subset by row and column.}
#'   \item{\code{subset}, \code{subsetProt} and \code{subsetPept}:}{Allows
#'   to subset based on a logical comparison using a column name from the
#'   rowData \code{data.frame}.}
#' }
#'
#' @section Merging:
#'
#' The following functions can be used to aggregate objects of the same class
#' together:
#'
#' \describe{
#'   \item{\code{cbind}:}{Joins two or more objects horizontally (adding
#'   samples). Must have the same proteins/peptides and in the same order.}
#'   \item{\code{rbind}:}{Joins two or more objects vertically (adding
#'   proteins/peptides). Must have the same samples and in the same order.}
#'   \item{\code{merge}:}{Joins two objects by adding new samples and tries
#'   to merge the proteins/peptide rowData data.frames and recalculate the
#'   linkerDf \code{data.frame} for the \code{ProteomicsExperiment} class.}
#' }
#'
#' @examples
#' # Accessors
#' ## assays
#' assays(wormsPE)
#' assaysProt(wormsPE)
#' assaysPept(wormsPE)
#'
#' ## assaysNames
#' assayNames(wormsPE)
#' assayNamesProt(wormsPE)
#' assayNamesPept(wormsPE)
#'
#' ## colData
#' colData(wormsPE)
#'
#' ## rowData
#' rowData(wormsPE)
#' rowDataProt(wormsPE)
#' rowDataPept(wormsPE)
#'
#' ## metadata
#' metadata(wormsPE)
#'
#' ## metaoptions
#' #metaoptions(wormsPE)
#'
#' ## linkerDf
#' linkerDf(wormsPE)
#'
#' # Dimensions
#' nrow(wormsPE)
#' ncol(wormsPE)
#' dim(wormsPE)
#' length(wormsPE)
#'
#' # Subsetting
#' wormsPE$line
#' #wormsPE[1,1]
#' #subsetProt(wormsPE, protein_id == 'AC3.2')
#' #subsetPept(wormsPE, Sequence == 'AIQEISDYHFLIK')
#'
#' # Merging
#' rbind(wormsPE[1:10, ], wormsPE[11:20, ])
#' cbind(wormsPE[,1:2], wormsPE[,3:4])
#' #merge(wormsPE[1:10, 1:3], wormsPE[3:10, 4:5])
NULL

