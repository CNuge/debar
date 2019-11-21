# debar
## A denoiser for DNA barcode data 
[![Build Status](https://travis-ci.com/CNuge/debar.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/debar)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
--------------------------------------------------------

`debar` is an R package designed to denoise COI-5P DNA barcode data generated through high throughput sequencing data. 


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

All a user needs to do is specify the input and output files and any custom paramaters (see `?denoise` or the manual for exhautive list). The `denoise_file` function accepts barcode data in either `fastq` or `fasta` formats (gzipped (`.gz`) files are also permitted).
Small example inputs are included with the package. After installation, these can be accessed as follows:
```
#fasta
fasta_example_file = system.file('extdata/coi_sequel_data_subset.fasta', package = 'debar')
#fastq
fastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq', package = 'debar')
#gzipped fastq
gzfastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq.gz', package = 'debar')
```

Note: running these examples will generate output files [in your current working directory!](https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces)
```
denoise_file(fastq_example_file, filename = "example_output.fastq")
```

File-to-file denoising can also be parallelized across multiple CPU cores. The denoising of each sequences in the input file is conducted separately, so using multiple cores will decrease the time needed to complete denoising roughly by a factor of the number of cores used. If you are denoising complete sequencer outputs, it is highly reccommend that you do so with as many cores as possible. For example, denoising of 160,000 sequence reads on a 64-core server (all default paramaters) takes approximately 1hr and 42 minutes, while on a laptop with only 8 cores would take almost 14 hours!

```
#debar works best when the tasks are highly parallelized
denoise_file(fastq_example_file, filename = "multicore-example_output.fastq", multicore = 8, log_file = TRUE, keep_rejects = TRUE) #set the multicore parameter to the number of CPU cores available
```

Certain paramater selections can further increase the speed with which `debar` can process data, but come with certain trade offs (that may or may not be worth consideration in the processing of your own data). The most drastic speed imrpovement is the provided 




### Denoising within R





## Version Notes

Initial design and default parameters are based on using seqDenoise with [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) on [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). Despite this, the package is designed to interface with fastq or fasta files of any origin (although the developers have yet to quantify performance on other data sources). The package uses a profile hidden Markov model (PHMMs) to identify and correct insertion and deletion errors within COI-5P sequences.


## Acknowledgements









