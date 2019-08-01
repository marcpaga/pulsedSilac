#' ProteomicsExperiment with pulsed silac data from C. elegans strains
#'
#' A pre-built \code{ProteomicsExperiment} object with data from a pulsed silac
#' experiment done in \emph{C. elegans} by Visscher et al. 2016. It only
#' contains the data from the first 250 priteins and two old worms strains
#' (OW40 and OW450).
#'
#' It is used as example in the pulsed silac vignette to illustrate the main
#' data analysis functions and in the examples of the documentation.
#'
#' @format A \code{ProteomicsExperiment} object with 250 proteins and 3574
#' peptides in a total of 14 samples.
#' \describe{
#'   \item{colData}{A \code{DataFrame} with the design of the experiment:
#'   samples, timepoints, replicates...}
#'   \item{assaysProt}{A list of matrices with quantification data at protein
#'   level: total intensity (int_total), light isotope intensity (int_light),
#'   heavy isotope intensity (int_heavy) and heavy/light isotope intensty ratio
#'   (ratio).}
#'   \item{rowDataProt}{A \code{DataFrame} with 22 columns that contains general
#'   protein information: ids, gene names, molecular weight...}
#'   \item{assaysPep}{A list of matrices with quantification data at peptide
#'   level: total intensity (int_total), light isotope intensity (int_light),
#'   heavy isotope intensity (int_heavy) and heavy/light isotope intensty ratio
#'   (ratio).}
#'   \item{rowDataPept}{A \code{DataFrame} with 46 columns that contains general
#'   protein information: ids, amino acids counts, length...}
#'   \item{linkerDf}{A \code{data.frame} with 3574 rows and 4 columns. It
#'   contains the relationships between proteins and peptides in the
#'   ProteomicsExperiment object.}
#' }
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/28679685}
'wormsPE'

#' ProteomicsExperiment with pulsed silac data from MEFs
#'
#' A pre-built \code{ProteinExperiment} object with data from a pulsed silac
#' experiment done in mouse embryonic fibroblasts (MEFs). Two cell cultures are
#' compared: cultured with or without serum.
#'
#' This dataset is used as an example, in the pulsed silac vignette, to show
#' the effect of comparing protein turnover between cell lines growing at
#' different rates.
#'
#' @format A \code{ProteinExperiment} object with 5223 proteins in a total of
#' 10 samples.
#' \describe{
#'   \item{colData}{A \code{DataFrame} with the design of the experiment:
#'   samples, timepoints, replicates...}
#'   \item{assays}{A list of matrices with quantification data at protein
#'   level: ratio and fraction.}
#'   \item{rowData}{A \code{DataFrame} with 3 columns that general protein
#'   id information.}
#' }
'mefPE'

#' Data frame with miss cleaved peptides quantifications from MaxQuant
#'
#' A \code{data.frame} that contains the output from searching heavy label
#' incorporation as amino acid modifications. This search was done using the
#' same data as in the \code{\link{wormsPE}} data.
#'
#' The reason why the search was done using isotopes as modifications is
#' because MaxQuant only looks for peptides in which all amino acids are
#' isotope labelled, miss-cleaved amino acids that contain an isotope mix are
#' not measured by the normal silac search engine. This allows to find peptides
#' which a mix of incoporated isotopes.
#'
#' This dataset is used as an example, in the pulsed silac vignette,
#' to estimate the amount of old isotope label in newly synthesized proteins
#' (amino acid recycling).
#'
#' @format A \code{data.frame}
#' \describe{
#'   \item{Sequence}{Column used as id to match with the peptides in the
#'   \code{\link{wormsPE}} ProteomicsExperiment.}
#'   \item{Modifications}{Column indicating how many heavy isotopes the peptide
#'   has.}
#'   \item{Intensity.*}{Columns containing quantification data in each sample.}
#'   \item{...}{Other columns containing peptide related information.}
#' }
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/28679685}
'recycleLightLysine'


