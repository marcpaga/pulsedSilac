# Documentation for the getters and setters for the ProteinExperiment,
# PeptideExperiment and ProteomicsExperiment classes

#' @rdname ProteinPeptideExperiment-accessors
#' @name ProteinPeptideExperiment-accessors
#'
#' @title Accessors for the ProteinExperiment and PeptideExperiment classes
#'
#' @description All the accessors, dimension, subsetting, merging and coercers
#' that work on \code{ProteinExperiment} and \code{PeptideExperiment} objects.
#' Functions that work on \code{SummarizedExperiment} objects should also work
#' on these two objects. Detailed examples of these functions can be found in
#' the vignette of this package.
#'
#' @param x A \code{ProteinExperiment} or a \code{PeptideExperiment} object.
#' @param ... For \code{rbind} and \code{cbind} are \code{ProteinExperiment} or
#' \code{PeptideExperiment} objects to be joined together. For \code{subset}
#' it is a logical comparison using a column name from the respective rowData
#' \code{data.frame}. Otherwise unused.
#' @param value An object of class specified in the S4 method signature or as
#' described in the following sections.
#' @param deparse.level Unused.
#'
#' @section Accessors:
#'
#' The following functions can be used to access the data in the class slots
#'
#' \describe{
#'   \item{\code{assays}:}{Access the assays (list of matrices) of the object.
#'   Value should be a matrix or list of matrices.}
#'   \item{\code{assayNames}}{ Access the assay names of the object. Value
#'   should be a character vector.}
#'   \item{\code{rowData}}{Access the protein/peptide feature \code{data.frame}
#'   of the object. Value should be \code{data.frame} with as many rows as
#'   proteins/peptides.}
#'   \item{\code{colData}:}{Access the samples \code{data.frame} of the object.
#'   Value should be a \code{data.frame} with as many rows as samples.}
#'   \item{\code{metadata}:}{Access the metadata \code{list} of the object.
#'   Value should be a \code{list}.}
#'   \item{\code{metaoptions}:}{Access the metaoptions \code{list} of the
#'   object. Value should be a \code{list}.}
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
#'   \item{\code{ncol}:}{Gives how many samples the object has.}
#'   \item{\code{dim}:}{Gives the number of proteins/peptides and the
#'   number of samples the object has.}
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
#'   \item{\code{subset}:}{Allows to subset based on a logical comparison
#'   using a column name from the rowData \code{data.frame}.}
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
#'   to merge the proteins/peptide rowData data.frames.}
#' }
#'
#' Merge methods are explained in detail in \link{merge}.
#'
#' @section Coercers:
#'
#' The folloing functions can be used to transform a \code{ProteinExperiment}
#' or a \code{PeptideExperiment} into a \code{SummarizedExperiment} or a
#' \code{data.frame}.
#'
#' \describe{
#'   \item{\code{as(x, 'SummarizedEpriment')}:}{Transforms the object into an
#'   object of class \code{SummarizedExperiment}.}
#'   \item{\code{as(x, 'data.frame')}:}{Transforms the object into an
#'   object of class \code{data.frame}.}
#' }
#'
#' @return Elements from \code{ProteinExperiment} or \code{PeptideExperiment}
#' objects.
#'
#' @examples
#' data('wormsPE')
#' protPE <- ProtExp(wormsPE)
#'
#' # Accessors
#' ## assays
#' assays(protPE)
#'
#' ## assaysNames
#' assayNames(protPE)
#'
#' ## colData
#' colData(protPE)
#'
#' ## rowData
#' rowData(protPE)
#'
#' ## metadata
#' metadata(protPE)
#'
#' ## metaoptions
#' #metaoptions(protPE)
#'
#'
#' # Dimensions
#' nrow(protPE)
#' ncol(protPE)
#' dim(protPE)
#' length(protPE)
#'
#' # Subsetting
#' protPE$line
#' protPE[1,1]
#' subset(protPE, protein_id == 'AC3.2')
#'
#' # Merging
#' rbind(protPE[1:10, ], protPE[11:20, ])
#' cbind(protPE[,1:2], protPE[,3:4])
#' #merge(protPE[1:10, 1:3], protPE[3:10, 4:5])
#'
#' # Coercers
#' as(protPE, 'SummarizedExperiment')
#' as(protPE, 'data.frame')
NULL



#' @rdname ProteomicsExperiment-accessors
#' @name ProteomicsExperiment-accessors
#'
#' @title Accessors for the ProteomicsExperiment class
#'
#' @description All the accessors, dimension, subsetting, merging and coercers
#' that work on \code{ProteomicsExperiment} objects. Since the
#' \code{ProteomicsExperiment} object has both protein and peptide level data,
#' most of the functions have a 'Prot' or 'Pept' suffix to indicate which
#' level should be used. If the non-suffix function is used, then a list with
#' both protein and peptide data is returned.
#' These functions also work on \code{ProteinExperiment}
#' and \code{PeptideExperiment} objects.
#'
#' @param x A \code{ProteomicsExperiment} object.
#' @param i,j For \code{`[`}, \code{i}, \code{j} are subscripts that can act to
#' subset the rows and columns of \code{x}.
#' @param drop A \code{logical} indicating if dimensions should be lowered if
#' possible when subsetting.
#' @param ... For \code{rbind} and \code{cbind} are \code{ProteinExperiment} or
#' \code{PeptideExperiment} objects to be joined together. For \code{subset}
#' it is a logical comparison using a column name from the respective rowData
#' \code{data.frame}. Otherwise unused.
#' @param value An object of class specified in the S4 method signature or as
#' described in the following sections.
#' @param name Column name of colData.
#' @param use.names Unused.
#' @param withDimnames Unused.
#' @param deparse.level Unused.
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
#'   \item{\code{ProtExp} and \code{PeptExp}:}{Access the experiment objects
#'   in a \code{ProteomicsExperiment}.}
#' }
#'
#' @section Dimensions:
#'
#' The following functions can be used to get the number of proteins/peptides
#' and number of samples:
#'
#' \describe{
#'   \item{\code{nrow}:}{Gives how many proteins and peptides the object
#'   has.}
#'   \item{\code{ncol}:}{Gives how many samples the object has.}
#'   \item{\code{dim}:}{Gives both the number of proteins andpeptides and the
#'   number of samples the object has.}
#'   \item{\code{length}:}{Gives how many proteins and or peptides the object
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
#' Merge methods are explained in detail in \link{merge}.
#'
#' @return Elements from a \code{ProteomicsExperiment} object.
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
#' wormsPE[1,1]
#' subsetProt(wormsPE, protein_id == 'AC3.2')
#' subsetPept(wormsPE, Sequence == 'AIQEISDYHFLIK')
#'
#' # Merging
#' rbind(wormsPE[1:10, ], wormsPE[11:20, ])
#' cbind(wormsPE[,1:2], wormsPE[,3:4])
#' #merge(wormsPE[1:10, 1:3], wormsPE[3:10, 4:5])
NULL



