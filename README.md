# pulsedSilac

## Summary

pulsedSilac is an R package developed to analyze mass spectrometry data using the SILAC (Stable Isotope Label by Amino acids in Cell culture) over time. Functions are provided to organize the data, calculate isotope ratios, isotope fractions, model protein turnover, compare turnover models, estimate cell growth and estimate isotope recycling. Several visualization tools are also included to do basic data exploration, quality control, condition comparison, individual model inspection and model comparison.

## Installation

The package can be installed from Bioconductor by running:

```{r}

if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("pulsedSilac")

```

To get the latest development version from Github, the devtools package can be used:

```{r}

if (!'devtools' %in% installed.packages()) {
  install.packages('devtools')
}
devtools::install_github("marcpaga/pulsedSilac")

```
