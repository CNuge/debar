
#' \pkg{debar}
#'
#' \pkg{debar} is an R package designed for the identification and removal 
#' of insertion and deletion errors from COI-5P DNA barcode data.
#' 
#' @details 
#' \pkg{debar} is built around the DNAseq object, which takes a COI-5P DNA barcode sequence
#' and optionally its associated name and PHRED quality information as input. The package utilizes
#' a nucleotide profile hidden Markov model (PHMM) for the identification of the COI-5P region of 
#' an input sequence and the identification and correction of indel errors from within the COI-5P 
#' region of the sequence. Indel corrections are by default applied in a conservative fashion, with subsequent
#' censorship of 7 base pairs in either direction of an indel correction to mask most instances where the exact
#' bp corresponding to the indel was not found exactly. Numerous filtering and double check steps are applied, and
#' the package includes functions for input/output for either fasta or fastq formats.
#' 
#' The denoise pipeline is heavily paramaterized so that a user can tailor the denoising execution for their own 
#' data structure and goal.
#' 
#' @section Functions:
#' \itemize{
#' \item \code{\link{denoise_file}} Run the denoise pipeline for each sequence in a specified input file.
#' \item \code{\link{denoise}} Run the denoise pipeline for a specified sequence
#' 
#' \item \code{\link{read_fasta}} Read data from a fasta file to a data frame.
#' \item \code{\link{read_fastq}} Read data from a fastq file to a data frame. 
#' \item \code{\link{write_fasta}} Write a denoised sequence to the specified fasta file.
#' \item \code{\link{write_fastq}} Write a denoised sequence and associated quality information to the specified fastq file.
#' 
#' \item \code{\link{DNAseq}} Builds a DNAseq class object
#' \item \code{\link{frame}} Match a sequence against the COI-5P PHMM using the 
#' Viterbi algorithm to establish the reading frame,
#' optional rejection of sequence based on the quality of the match to the PHMM.
#' \item \code{\link{adjust}} Use the PHMM path output corresponding to the sequence to 
#' adjust the DNA sequence and remove indels.
#' Optional censorship of sequence around the corrections.
#' \item \code{\link{aa_check}} Translate the adjusted sequence to amino acids and check it for stop codons.
#' \item \code{\link{outseq}} Construct the output data for the given sequence. 
#' Optionally can include or exculde sequence data from
#' outside of the COI-5P region (part of sequence that was not denoised).
#' }
#' 
#' @section Data:
#' \itemize{
#' \item \code{\link{example_nt_string}} An example COI-5P sequence with no errors.
#' \item \code{\link{example_nt_string_errors}} An example COI-5P sequence with two indel errors.
#' }
#' @author Cameron M. Nugent
#' @docType package
#' @name debar
#################
NULL
