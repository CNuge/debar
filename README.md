# debar
## A denoiser for DNA barcode data 
[![Build Status](https://travis-ci.com/CNuge/debar.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/debar)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
--------------------------------------------------------

`debar` is an R package designed to denoise COI-5P DNA barcode data generated through high throughput sequencing data. 


## Installation




## Use and examples

The package's vignette contains detailed explinations of the functions and parameters of the `debar` denoising pipeline. The use is encouraged to read this document in order to get oriented and effectively deploy `debar` in the denoising of their own data. Following package installation, the vignette can be accessed from within R through the following command:
```
vignette('debar-vignette')
```


### File to file denoising


### Denoising within R





## Version Notes

Initial design and default parameters are based on using seqDenoise with [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) on [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). Despite this, the package is designed to interface with fastq or fasta files of any origin (although the developers have yet to quantify performance on other data sources). The package uses a profile hidden Markov model (PHMMs) to identify and correct insertion and deletion errors within COI-5P sequences.


## Acknowldgements









