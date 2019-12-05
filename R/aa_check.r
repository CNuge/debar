
#' Translate the sequence and it for stop codons
#'
#' A side product of the framing and adjustment functions is that the reading frame of the sequence is established
#' and translation can be conducted with high confidence. If the adjustment did in fact correct indel errors, the
#' translated sequence should feature no stop codons. If stop codons are present this is grounds for rejection as 
#' indel errors (or some other large issue) persists in the sequence.
#' 
#' This test has limitations, as any indels late in the DNA sequence may not lead to stop codons existing. Additionally
#' by default censored translation is used by this function when producing the amino acid sequence, so as to 
#' eliminate taxonomic bias against organisms with esoteric translation tables. The likelihood of catching errors is
#' increased if the genetic code corresponding to sequences is known.
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
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @examples
#' #previously called
#' ex_data = DNAseq(example_nt_string, name = 'ex1')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data)
#' #run the aa check on the adjusted sequences
#' ex_data = aa_check(ex_data)
#' ex_data$aaSeq #view the amino acid sequence
#' ex_data$stop_codons # Boolean indicating if stop codons in the amino acid sequence
#' 
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
    #turn the dashes into n so that we can interface with sequinr
    dna_list = gsub('-', 'n', as.character(tolower(x$adjusted_sequence)))
    #translate using the designated numcode, returns a vector of AAs
    aa_vec = seqinr::translate(dna_list, frame = 0, numcode=trans_table, ambiguous= TRUE, NAstring = '-')
    x$aaSeq = paste(aa_vec, collapse= "")
  }
  
  #if there is a stop codon
  if(grepl('\\*', x$aaSeq)){
    x$stop_codons = TRUE
    x$reject = TRUE
    # below is for if the sequence should be masked in case of stop codons
    # this functionality currently disabled - not conservative enough
    #censor_at = first_stop(x$aaSeq)
    #censor_str = strsplit(x$adjusted_seq, "")[[1]]
    #censor_str[censor_at:length(censor_str)] = "N"
    #x$adjusted_seq = paste(censor_str, collapse = "")
  }else{
    x$stop_codons = FALSE
  }
  
  return(x)  
}

