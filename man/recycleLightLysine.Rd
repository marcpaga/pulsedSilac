% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{recycleLightLysine}
\alias{recycleLightLysine}
\title{Data frame with miss cleaved peptides quantifications from MaxQuant}
\format{A \code{data.frame}
\describe{
  \item{Sequence}{Column used as id to match with the peptides in the
  \code{\link{wormsPE}} ProteomicsExperiment.}
  \item{Modifications}{Column indicating how many heavy isotopes the peptide
  has.}
  \item{Intensity.*}{Columns containing quantification data in each sample.}
  \item{...}{Other columns containing peptide related information.}
}}
\usage{
data(recycleLightLysine)
}
\description{
A \code{data.frame} that contains the output from searching heavy label
incorporation as amino acid modifications. This search was done using the
same data as in the \code{\link{wormsPE}} data.
}
\details{
The reason why the search was done using isotopes as modifications is
because MaxQuant only looks for peptides in which all amino acids are
isotope labelled, miss-cleaved amino acids that contain an isotope mix are
not measured by the normal silac search engine. This allows to find peptides
which a mix of incoporated isotopes.

This dataset is used as an example, in the pulsed silac vignette,
to estimate the amount of old isotope label in newly synthesized proteins
(amino acid recycling).
}
\references{
\url{https://www.ncbi.nlm.nih.gov/pubmed/28679685}
}
\keyword{datasets}
