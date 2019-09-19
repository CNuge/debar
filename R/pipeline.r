#' Run the denoiser pipeline for a group of ccs reads.
#'
#' 
#' @param x a ccs_reads class object.
#' @param ... additional arguments to be passed between methods.
#' @param censor_length the number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 3.
#' @param ambig_char The character to use for ambigious positions in the sequence.
#' @param to_file Boolean indicating whether the consensus sequences should be written to a file. Default is TRUE.
#' @param outformat The format of the output file. Options are fasta (default) or fastq format.
#' @param filename The name of the file to output the data to. Default filenames are: denoised.fasta or denoised.fastq.
#' @param append Should the ccs consensus sequence be appended to the output file?(TRUE) 
#' Or overwrite the file?(FALSE) Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'. Used with write_fastq only
#' @param aa_check Boolean indicating whether the amino acid sequence should be generated and assessed for each ccs read. Default = FALSE.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids. Default is "auto", meaning
#' that the translation table will be inferred from the ccs_reads object's order. Used only when aa_check = TRUE.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would offset the sequence by one and therefore make the second character in the
#' framed sequence the the first nucelotide of the first codon. Used only when aa_check = TRUE.
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{build_ccs}}
#' @examples
#' # construct the ccs_reads object
#' ex_data = build_ccs(ex_ccs_read_list, order = 'Diptera', id = 'SSGBC787-14')
#' # denoise the ccs_reads object
#' ex_data = denoise(ex_data)
#' @export
#' @name denoise
denoise = function(x, ...){
  UseMethod("denoise")
}


#' @rdname denoise
#' @export
denoise.ccs_reads = function(x, ...,
                             ambig_char = "N",
                             censor_length = 3,
                             to_file = TRUE,
                             outformat = "fasta", 
                             filename = NULL, 
                             phred_placeholder = "#",
                             aa_check = FALSE, 
                             trans_table = "auto",
                             frame_offset = 0,
                             append = TRUE
                             ){
  if(outformat != "fastq" && outformat != "fasta" && outformat != "none"){
    stop("Invalid output format! Must be one of: 'fasta', 'fastq' or 'none'")
  }
  
  x = frame(x)
  x = adjust(x, censor_length = censor_length)
  if(aa_check == TRUE){
    x = aa_check(x, trans_table = trans_table, frame_offset = frame_offset)
  }
  x = consensus(x, ambig_char = ambig_char)
 
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
        write_fastq(x, phred_placeholder = phred_placeholder, append = append)
      }else{
        write_fastq(x, filename = filename, phred_placeholder = phred_placeholder, append = append)
      }
    }    
  } 
  x
}

