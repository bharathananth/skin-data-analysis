<!-- badges: start -->
<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205155"><img src="https://img.shields.io/badge/Data-GSE205155-green.svg?style=plastic" alt="" /></a>
[![](https://img.shields.io/badge/doi-10.1101/2022.06.03.494693-yellow.svg)](https://doi.org/10.1101/2022.06.03.494693)
[![DOI](https://zenodo.org/badge/541140885.svg)](https://zenodo.org/badge/latestdoi/541140885)
 <!-- badges: end -->

# Skin-data-analysis


This repository contains the reproducible code for the analysis and figures in the manuscript "Inter-layer and inter-subject variability of circadian gene expression in human skin" ([preprint](https://doi.org/10.1101/2022.06.03.494693)) . The dependencies of the code are handled using the package `renv`.

To execute this code:


1. Clone this project at a suitable location (this places the code within a directory **skin-data-analysis**). 
2. Install the package `renv` with `install.packages("renv")`.
3. Open project `skin-data-analysis.Rproj` within directory **skin-data-analysis**.
4. Activate `renv` using `renv::activate()`.
5. Install the correct package versions into the local R library using `renv::restore()` (this can sometimes take quite a long time).

The R files can now be executed in order *0_..*, *1_..*, ... within the project.
