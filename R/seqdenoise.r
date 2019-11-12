
#' \pkg{seqdenoise}
#'
#' \pkg{seqdenoise} is an R package designed for the identification and removal 
#' of insertion and deletion errors from COI-5P DNA barcode data.
#' 
#' @details 
#' \pkg{seqdenoise} is built around the DNAseq object, which takes a COI-5P DNA barcode sequence
#' and optionally its associated name and PHRED quality information as input. The package utilizes
#' a nucleotide profile hidden Markov model (PHMM) for the identification of the COI-5P region of 
#' an input sequence and the identification and correction of indel errors from within the COI-5P 
#' region of the sequence. Indel corrections are by default applied in a conservative fashion, with subsequent
#' censorship of 7 base pairs in either direction of an indel correction to mask most instances where the exact
#' bp corresponding to the indel was not found exactly. Numerous filtering and double check steps are applied, and
#' the package includes functions for input/output for either fasta or fastq formats.
#' 
#' The denoise pipeline is heavily paramaterized so that a user can tailor the denoising execution for their own 
#' data structure and gold.
#################
NULL
