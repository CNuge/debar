#' Run the denoiser pipeline for a sequence read.
#'
#' 
#' @param x a DNA sequence string.
#' @param ... additional arguments to be passed between methods.
#' @param name an optional character string. Identifier for the sequence.
#' @param censor_length the number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 5.
#' @param to_file Boolean indicating whether the sequence should be written to a file. Default is TRUE.
#' @param ambig_char The character to use for ambigious positions in the sequence that is output to the file. Default is N.
#' @param keep_flanks Should the regions of the input sequence outside of the barcode region be readded to the denoised sequence
#' prior to outputting to the file. Default is TRUE. 
#' False will lead to only the denoised sequence for the 657bp barcode region being output to the file.
#' @param outformat The format of the output file. Options are fasta or fastq (default) format.
#' @param filename The name of the file to output the data to. Default filenames are respectively: denoised.fasta or denoised.fastq.
#' @param append Should the denoised sequence be appended to the output file?(TRUE) 
#' Or should the sequence overwrite the output file?(FALSE) Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'. Used with write_fastq only.
#' @param aa_check Boolean indicating whether the amino acid sequence should be generated and assessed for stop codons. Default = FALSE.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids. Default is 0, meaning
#' that censored translation is performed (amigious codons ignored). Used only when aa_check = TRUE.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would offset the sequence by one and therefore make the second character in the
#' framed sequence the the first nucelotide of the first codon. Used only when aa_check = TRUE.
#' @return a class object of code{"DNAseq"} 
#' @examples
#' # Denoise example sequence with default paramaters.
#' ex_data = denoise(example_nt_string, name = 'example_sequence_1')
#' @export
#' @name denoise
denoise = function(x, ...){
  UseMethod("denoise")
}


#' @rdname denoise
#' @export
denoise = function(x, ...,
                             name = character(),
                             ambig_char = "N",
                             censor_length = 5,
                             to_file = TRUE,
                             keep_flanks = TRUE,
                             outformat = "fastq", 
                             filename = NULL, 
                             phred_placeholder = "#",
                             aa_check = FALSE, 
                             trans_table = 0,
                             frame_offset = 0,
                             append = TRUE
                             ){
  if(outformat != "fastq" && outformat != "fasta" && outformat != "none"){
    stop("Invalid output format! Must be one of: 'fasta', 'fastq' or 'none'")
  }
  
  x = DNAseq(x, name = name )
  x = frame(x)
  x = adjust(x, censor_length = censor_length)
  if(aa_check == TRUE){
    x = aa_check(x, trans_table = trans_table, frame_offset = frame_offset)
  }
  
  if(to_file == TRUE){
    if(outformat == "fasta"){
      if(is.null(filename)){
        write_fasta(x, append = append)
      }else{
        write_fasta(x, filename = filename, append = append)
      }
    }
    
    if(outformat == "fastq"){
      if(is.null(filename)){
        write_fastq(x, ambig_char= ambig_char,
                       keep_flanks = keep_flanks, 
                       phred_placeholder = phred_placeholder, 
                       append = append)
      }else{
        write_fastq(x, ambig_char = ambig_char,
                       keep_flanks = keep_flanks, 
                       filename = filename, 
                       phred_placeholder = phred_placeholder, 
                       append = append)
      }
    }    
  } 
  x
}


#'Denoise the sequence data from a given file
#'@param filename The name of the file to
#'@param file_type
#'@param multicore
#'@param ... additional arguments to be passed to the denoise function.
#'
#'@seealso \code{\link{denoise}}
#'  
denoise_file = function(filename, file_type = "fastq", multicore = FALSE, ...){
  if(file_type == "fastq"){
    data = read_fastq(filename)
  }else if (file_type == "fasta"){
    sequences = read_fasta(filename)
  }else{
    stop("file_type must be either fasta or fastq")
  }
  
  if(multicore == FALSE){
    for(i in 1:length(sequences)){
      denoise(data$sequence[[i]], name = data$header_data[[i]], ...)
    }
  }
}
