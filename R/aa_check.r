#TODO - possibly shift the aa_check application down to the level of adjusted sequences
# and have it be applied prior to the consensus step. Probably more useful in that capacity.


#' build an AAbin with ape.
#'
#' @keywords internal
individual_AAbin = function(aa_string){
  return(ape::as.AAbin(strsplit(as.character(aa_string),"")))
}

#' Find the first stop codon in the amino acid sequence. 
#' Return the corresponding nucleotide mask start position.
#'
#' @keywords internal
first_stop = function(aa_str){
  (nchar(strsplit(aa_str, '\\*')[[1]][[1]])+1)*3
}

#' Translate the sequence and check log likelihood as a final double check.
#'
#' 
#' @param x a ccs_reads class object.
#' @param ... additional arguments to be passed between methods.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids. Default is "auto", meaning
#' that the translation table will be inferred from the ccs_reads object's order.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would offset the sequence by one and therefore make the second character in the
#' framed sequence the the first nucelotide of the first codon.
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{build_ccs}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @examples
#' #previously called
#' ex_data = build_ccs(ex_ccs_read_list, order = 'Diptera', id = 'SSGBC787-14')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data)
#' ex_data = consensus(ex_data)
#' #run the aa check on the adjusted sequences
#' ex_data = aa_check(ex_data)
#' @export
#' @name aa_check
aa_check = function(x, ...){
  UseMethod("aa_check")
} 

#' @rdname aa_check
#' @export
aa_check.DNAseq = function(x, ..., trans_table = 0, frame_offset = 0){
  
  if(trans_table == 0){
    x$aaSeq = censored_translation(paste(x$adjusted_sequence, collapse=""), reading_frame = (frame_offset+1))
  }else{
    #split the DNA string into a vector, all characters to lower case
    dna_list = strsplit(gsub('-', 'n', as.character(tolower(x$adjusted_sequence))),"")
    dna_vec = dna_list[[1]]
    #translate using the designated numcode, returns a vector of AAs
    aa_vec = seqinr::translate(dna_vec, frame = frame_offset, numcode=trans_table, ambiguous= TRUE, NAstring = '-')
    
    x$aaSeq = paste(aa_vec, collapse= "")
  }
  
  #if there is a stop codon
  if(grepl('\\*', x$aaSeq)){
    x$stop_codons = TRUE
    censor_at = first_stop(x$aaSeq)
    censor_str = strsplit(x$adjusted_seq, "")[[1]]
    censor_str[censor_at:length(censor_str)] = "N"
    x$adjusted_seq = paste(censor_str, collapse = "")
  }else{
    x$stop_codons = FALSE
  }
  
  return(x)  
}

