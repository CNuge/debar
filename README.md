# debar
## DEnoise BARcode data with profile hidden Markov models 
[![Build Status](https://travis-ci.com/CNuge/debar.svg?token=H6eQaqsE1kLqYX3zZ1Xz&branch=master)](https://travis-ci.com/CNuge/debar)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![codecov](https://codecov.io/gh/CNuge/debar/branch/master/graph/badge.svg)](https://codecov.io/gh/CNuge/debar)
--------------------------------------------------------

`debar` is an R package designed for denoising sequence data for the animal DNA barcode marker: cytochrome c oxidase I (COI-5P, or the five prime portion of COI). The package is designed to detect and correct insertion and deletion errors within barcode sequences. This is accomplished through comparison of input sequences against a profile hidden Markov (PHMM) model [using the Viterbi algorithm](https://en.wikipedia.org/wiki/Viterbi_algorithm) and adjustment of the sequence based on the reported Viterbi path (`debar` depends on functions from the R package [aphid](https://CRAN.R-project.org/package=aphid) for the PHMM strucutre and for running the Viterbi algorithm).

Inserted base pairs are removed and deleted base pairs are accounted for through the introduction of a placeholder character. Since the PHMM is a probabilistic representation of the COI barcode, corrections are not always perfect. For this reason `debar` censors base pairs adjacent to reported indel sites, turning them into placeholder characters (default is 7bp in either direction, this feature can be disabled). Testing has shown that this censorship results in the correct sequence length being restored and erroneous base pairs being masked the vast majority of the time (>95%). Multiple denoised and censored reads from a sample can then be combined to obtain denoised barcode. 

## Intended use

The `debar` denoising method is designed to increase the accuracy of reported barcode sequences. It is therefore best applied when the accuracy of a barcode is paramount, namely in the generation of novel barcode sequences. In the case of metabarcoding analysis `debar` can allow for denoised halpotypes to be obtained for operational taxonomic units (OTUs).

The package was initially designed for processing [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) outputs produced by [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). However, the method is barcode-based (as opposed to sequencer-based) and should therefore be robust to COI sequence data of any origin. If you intend to apply debar in the denoising of metabarcode data, it is recommended that you do so after quality filtering and dereplication of reads. Please consult the package vignette for recommended parameters and examples of integration into barcoding and metabarcoding workflows. 

## Installation

The development version of `debar` can be installed directly from GitHub. You'll need to have the R package `devtools` installed and loaded. Also note if the build_vignettes option is set to true, you will need to have the R package `knitr` installed.

```
#install.packages("devtools")
#install.packages("knitr") #required if build_vignettes = TRUE
#library(devtools) 
devtools::install_github("CNuge/debar", build_vignettes = TRUE)
library(debar)
```
## Use and examples

The package's vignette contains detailed explanations of the functions and parameters of the `debar` denoising pipeline. The use is encouraged to read this document in order to get oriented and effectively deploy `debar` in the denoising of their own data. Following package installation, the vignette can be accessed from within R through the following command:
```
vignette('debar-vignette')
```
A second vignette, with a detailed explination of the denoising process is also included:
```
vignette('debar-algorithm-details')
```

### Denoising within R

`debar` can be used to perform denoising of sequences from within R. The `denoise` function can be used to process a single read and (optionally) its associated quality information as well.
```
#read a file of example sequences 
#fastq
fastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq', package = 'debar')
data = read_fastq(fastq_example_file)

#denoise a given read
denoised_seq = denoise(data$sequence[[1]], 
                      name = data$header_data[[1]],
                      quality = data$quality[[1]], 
                      to_file = FALSE)
?denoise # for an exhaustive list of parameter options.
names(denoised_seq) # for list of available object components

```
This will produce a DNAseq object, from which detailed information related to a given read can be accessed using the dollar sign notation. 

#### Batch processing
A list of sequences for a given haplotype can be denoised using the `denoise_list` function. This is an especially useful feature for tasks such as denoising sequences in an OTU prior to determining the consensus sequence, or in obtaining common haplotypes for an OTU. The `denoise_list` function is parallelized and can be run acrosss multiple cores (specify the number available with the `cores` argument)

```
# ex_nt_list is an example list of four barcode sequences that contain errors.
ex_out = denoise_list(ex_nt_list, cores = 2)
```

Optionally, when multiple sequences are available from a given sample (or OTU) they can be denoised as a group. Passing the `keep_flanks=FALSE` option to the function will produce denoised outputs with a common reading frame (leading placeholder Ns are added as needed). The `consensus_sequence` function can then be used to obtain a consensus from the denoised sequences.
```
ex_out = denoise_list(ex_nt_list, keep_flanks=FALSE)
ex_out #each output individually has some missing information
barcode_seq = consensus_sequence(ex_out)
barcode_seq #aligned through the denoising process, a consensus without missing information can be obtained
```


### File-to-file denoising
Denoising of COI-5P barcode data with `debar` can be conducted in a file-to-file fashion using the `denoise_file` function. 

All a user needs to do is specify the input and output files and as well as any deviations from the default parameters they wish to apply (see `?denoise` or the manual for exhaustive parameter list). The sequences in the input file will be denoised and written to the output file in the specified format. The `denoise_file` function accepts barcode data in either `fastq` or `fasta` formats (gzipped (`.gz`) files are also permitted). Small example inputs are included with the package. 

A complete file can be denoised in a single line of R code, simply specify the input and output files. When processing fastq files, `debar` preserves the phred scores.

*Note*: running the following example will generate an output file [in your current working directory!](https://support.rstudio.com/hc/en-us/articles/200711843-Working-Directories-and-Workspaces)
```
#gzipped fastq
gzfastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq.gz', package = 'debar')

denoise_file(gzfastq_example_file, outfile = "example_output.fastq")
```
If you are planning on utilizing `debar` for large input files, please consult the package vignette section 'Parameter combinations - speed and accuracy trade-offs' for suggestions on how to optimize performance when scaling to tens or hundreds of thousands of sequences.

## Version Notes

Initial design and default parameters are based on using `debar` to process the circular consensus sequences of [single molecule real-time (SMRT) sequencing](https://www.pacb.com/smrt-science/smrt-sequencing/) outputs produced by [the Pacific Biosciences SEQUEL platform](https://www.pacb.com/products-and-services/sequel-system/). Despite this, the package is designed to interface with fastq or fasta files of any origin (although the developers have yet to quantify performance on other data sources).
In the future we hope to quantify performance and provide informed hyper-parameter choices for outputs from other sequencing platforms. If you are interested in beta testing `debar` on barcode or metabarcode data from other platforms, please [contact Cam](https://cnuge.github.io), we would be happy to work with you to optimize `debar`'s functionality for other sequencing platforms.

## Acknowledgements

Funding for the development of this software was provided by grants in Bioinformatics and Computational Biology from the Government of Canada through Genome Canada and Ontario Genomics and from the Ontario Research Fund. Funders played no role in the study design or preparation of this software.