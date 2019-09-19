
#' Determining the most common character in each column.
#'
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



#' Take the list of adjusted sequences and obtain the consensus sequence.
#'
#' @param x The list of adjusted DNA sequences, post censorship and AA correction
#'
#' @keywords internal
consensus_sequence = function(x){
  #pad the sequences so that they're all the same length
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
  seq_mat[seq_mat == "-"] = NA
  
  consensus_vec = apply(seq_mat, 2, col_mode)
  
  return(toupper(paste(consensus_vec, collapse="")))
}



#' Get the consensus sequence of the denoised circular consensus reads.
#'
#' This function takes the adjusted sequence reads and determines the consensus sequence.
#' Placeholder 'N' is added for instances where the nucleotide for the given sequence position 
#' cannot be stated with certainty. 
#' @param x a ccs_reads class object.
#' @param ... additional arguments to be passed between methods.
#' @param ambig_char The character to use for ambigious positions in the sequence.
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{build_ccs}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @examples
#' #previously run
#' ex_data = build_ccs(ex_ccs_read_list, order = 'Diptera', id = 'SSGBC787-14')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data, censor_length = 3)
#' #get the consensus of the adjusted sequences.
#' ex_data = consensus(ex_data, ambig_char = 'N')
#' @export
#' @name consensus
consensus = function(x, ...){
  UseMethod("consensus")
} 

#' @rdname consensus
#' @export
consensus.ccs_reads = function(x, ..., ambig_char = "N"){
  
  #take the stacked adjusted reads and turn them into a matrix.
  #make all the dashes into NAs and then take the majority of
  #each column, excluding the NAs, if there are only NAs for a 
  #column then put the ambig_char in that position.
  
  #then paste the sequence into a single string and save it to x$denoised
  x$consensus = consensus_sequence(x$adjusted_sequences)
  if(ambig_char != "N"){
    x$consensus = gsub("N", ambig_char, x$consensus )
  }
  return(x)
}

