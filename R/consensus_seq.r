#' Determining the most common character in each column.
#'
#' @return A single character string.
#' @keywords internal
col_mode = function(col_data) {
  if (!is.factor(col_data)){
    col_data <- factor(col_data)
  } 
  
  tab_col = tabulate(col_data)
  
  if(max(tab_col) == 0){
    return("n")
  }
  
  consensus = levels(col_data)[tab_col == max(tab_col)]
  
  if (length(consensus)>1) {
    return("n")
  } 
  
  return(consensus)
}



#' Take the list of denoised sequences and obtain the consensus sequence.
#'
#' @param x The list of adjusted DNA sequences, post censorship and AA correction
#' @seealso \code{\link{denoise_list}}
#' @examples
#' #denoise list of sequences with the k
#' ex_out = denoise_list(ex_nt_list, keep_flanks=FALSE)
#' 
#' barcode_seq = consensus_sequence(ex_out)
#' @return A string containing the consensus nucleotide sequence.
#' @export
#' @name consensus
consensus_sequence = function(x){
  #pad the sequences so that they're all the same length
  x = lapply(x, function(i) strsplit(i, "")[[1]])
  max_len = 0
  for(seq in x){
    new_len = length(seq)
    if(max_len < new_len){
      max_len = new_len
    }
  }

  for(i in 1:length(x)){
    length(x[[i]]) = max_len

  }

  seq_mat  = do.call(rbind, x)
  seq_mat[seq_mat == "N"] = NA
  
  consensus_vec = apply(seq_mat, 2, col_mode)
  
  return(toupper(paste(consensus_vec, collapse="")))
}
