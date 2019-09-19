# coiDenoiser
## An R package for the denoising coi5p barcode data generatd through Pacific Biosciences single molecule, real-time (SMRT) sequencing.
[![Build Status](https://travis-ci.com/CNuge/coiDenoiser.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/coiDenoiser)
--------------------------------------------------------

`coiDenoiser` is an R package designed to denoise DNA barcode data generated through [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) on [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). The package uses taxonomic specific profile hidden Markov models (PHMMs) to identify and correct insertion and deletion errors within individual circular consensus sequences (CCS). When multiple CCSs are avaliable for a sample, the corrected sequences are then evaluated to determine the most likely true sequence.

## Installation
You can download the development version of `coiDenoiser` directly from GitHub. The process requires having the R package `devtools` installed and loaded.
```
#install.packages("devtools")
#library(devtools)
devtools::install_github("CNuge/coiDenoiser", build_vignettes = TRUE)
library(coiDenoiser)
```
The vignette demonstrating how to use coiDenoiser can be accessed using the following R command:
```
vignette("coiDenoiser-vignette")
```

## Usage





