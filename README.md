# debar
## DEnoise BARcode data with profile hidden Markov models 
[![Build Status](https://travis-ci.com/CNuge/debar.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/debar)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
--------------------------------------------------------

`debar` is an R package designed for denoising high throughput sequencing data for the animal DNA barcode marker cytochrome c oxidase I (COI-5P, or the five prime portion of COI). It is designed to detect and correct insertion and deletion errors within sequencer outputs. This is accomplished through comparison of input sequences against a profile hidden Markov model [using the Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm). The Viterbi path produces is used to identify where insertion or deletion errors exist. Without going into too much deail about the inner working of the package, this is how the package corrects indels:

```
#Correcting insertions
#the 'G' at position 6 is an insertion
"ATGATGATG"
#the sequence is compared to the PHMM, 
#the 'G' is assigned a hidden state of 'insert', indicated by a 2 in the path
"111112111"
#debar then alters the sequence, removing the G
"ATGATATG"

#Correcting deletions
#There is a missing base pair after the first T in the sequence
"AAAATGATGATG"
"111101111111"
AAAAT-GATGATG-"
```

## Installation

The development version of `debar` can be installed directly from GitHub. You'll need to have the R package `devtools` installed and loaded. Also note if the build_vignettes option is set to true, you will need to have the R package `knitr` installed.

```
#install.packages("devtools")
#install.packages("knitr") #required if build_vignettes = TRUE
#library(devtools) 
devtools::install_github("CNuge/debar", build_vignettes = TRUE)
library(coil)
```
## Use and examples

The package's vignette contains detailed explanations of the functions and parameters of the `debar` denoising pipeline. The use is encouraged to read this document in order to get oriented and effectively deploy `debar` in the denoising of their own data. Following package installation, the vignette can be accessed from within R through the following command:
```
vignette('debar-vignette')
```

### File-to-file denoising
Denoising of COI-5P data with `debar` can be conducted in a file-to-file fashion using the `denoise_file` function. 

All a user needs to do is specify the input and output files and specify any deviations from the default paramaters (see `?denoise` or the manual for exhautive list) and all the sequences in the input file will be denoised and written to the output file. The `denoise_file` function accepts barcode data in either `fastq` or `fasta` formats (gzipped (`.gz`) files are also permitted). Small example inputs are included with the package. After installation, these example files can be accessed as follows:

```
#fasta
fasta_example_file = system.file('extdata/coi_sequel_data_subset.fasta', package = 'debar')
#fastq
fastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq', package = 'debar')
#gzipped fastq
gzfastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq.gz', package = 'debar')
```
Note: running the following examples will generate output files [in your current working directory!](https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces)

A complete file can be denoised in a single line of R code, simply specify the input and output files.
```
denoise_file(fastq_example_file, filename = "example_output.fastq")
```
If you are planning on utilizing `debar` on compete sequencing runs, please consut the package vignette section 'Paramater combinations - speed and accuracy tradeoffs' for suggestions on how to optimize performance in different situations.

### Denoising within R

The file-to-file denoising method should serve the purposes of the majority of users. However, you are also able to perform denoising of sequences from within R (for the purposes of parameter tuning, extraction of additional data etc.).


```
#read file
```

```
#denoise
?denoise # for an exhaustive list of parameter options.
```

The individual components of the denoise function can each be accessed and called individually. An detailed explination of the denoising steps is provided in the package's vignette.

## Version Notes

Initial design and default parameters are based on using `debar` to process [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) outputs produced by [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). Despite this, the package is designed to interface with fastq or fasta files of any origin (although the developers have yet to quantify performance on other data sources). The package uses a profile hidden Markov model (PHMMs) to identify and correct insertion and deletion errors within COI-5P sequences. In the future we hope to quantify performance and provide informed hyperparamater choices for outputs from other sequencing platforms. If you are interested in beta testing `debar` on barcode or metabarcode data from other platforms, please [contact Cam](https://cnuge.github.io), we would be happy to work with you to optimize `debar`'s functionality for other sequencing platforms.

## Acknowledgements

Funding for the development of this software was provided by grants in Bioinformatics and Computational Biology from the Government of Canada through Genome Canada and Ontario Genomics and from the Ontario Research Fund. Funders played no role in the study design or preparation of this software.
