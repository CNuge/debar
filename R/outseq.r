
#' Get the final denoised output sequence for a read.
#'
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param keep_flanks Default is TRUE.
#' @param ambig_char The character to use for ambigious positions in the sequence.
#' @param adjust_limit the maximum number of corrections that can be applied to a sequence read. If this number is exceeded 
#' then the entire read is masked with ambigious characters. Default is 5.
#'
outseq = function(x, ...){
  UseMethod("outseq")
}

#' @rdname outseq
#' @export
outseq.DNAseq = function(x, ..., keep_flanks = TRUE, ambig_char = "N", adjust_limit = 5){
  
  if(keep_flanks == TRUE){
    
    dashes_rm = sum(c(length(x$frame_dat$front), as.integer(x$data$len_first_front), 1), na.rm = TRUE)
    if(dashes_rm > length(x$adjusted_sequence)){
      x$adjusted_sequence = NULL
    }
    
    #here the good part of the read is pasted back together with the flanking information
    #that was omitted.
    x$outseq = paste(c(x$data$raw_removed_front, #part removed in front of first frame() call
                       x$frame_dat$removed_lead, #part removed in second frame() call
                       x$adjusted_sequence[dashes_rm:length(x$adjusted_sequence)], #the 'sweet spot' that is denoised
                       x$data$adjusted_trimmed, # part removed in adj_seq
                       x$frame_dat$removed_end, #part removed in second frame() call
                       x$data$raw_removed_end), # part removed from the end of first frame() call
                     collapse = "")
    
    if(!is.null(names(x$adjusted_sequence)) || !is.null(names(x$frame_dat$removed_lead))|| !is.null(names(x$data$raw_removed_front))){
      x$outphred = paste(c(names(x$data$raw_removed_front), #part removed in front of first frame() call
                           names(x$frame_dat$removed_lead), #part removed in second frame() call
                           names(x$adjusted_sequence[dashes_rm:length(x$adjusted_sequence)]), #the 'sweet spot' that is denoised
                           names(x$data$adjusted_trimmed), # part removed in adj_seq
                           names(x$frame_dat$removed_end), #part removed in second frame() call
                           names(x$data$raw_removed_end)), # part removed from the end of first frame() call
                         collapse = "")
    }
  }else{
    x$outseq = paste(x$adjusted_sequence, collapse = "")
    if(!is.null(names(x$adjusted_sequence))){
      x$outphred = names(x$adjusted_sequence)
    }
  }
  
  x$outseq = toupper(gsub("-", ambig_char, x$outseq )) 
  
  if(x$adjustment_count > adjust_limit){
    x$reject = TRUE
  }else{
    x$reject = FALSE
  }
  
  return(x)
}