
#' Get the final denoised output sequence for a read.
#'
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param keep_flanks Should the regions of the input sequence outside of the barcode region be readded to the denoised sequence
#' prior to outputting to the file. Options are TRUE, FALSE and 'right'. The 'right' option will keep the trailing flank
#' but remove the leading flank. Default is TRUE. 
#' False will lead to only the denoised sequence for the 657bp barcode region being output to the file.
#' @param ambig_char The character to use for ambigious positions in the sequence.
#' @param adjust_limit the maximum number of corrections that can be applied to a sequence read. If this number is exceeded 
#' then the entire read is masked with ambigious characters. Default is 5.
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @examples
#' #previously run
#' excess_string = paste0("CCCCCC", example_nt_string_errors, 
#'                        "CCCCCCCC", collapse="")
#' ex_data = DNAseq(excess_string, name = 'ex1')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data)
#' #build output sequence with trimmed edges
#' ex_data = outseq(ex_data, keep_flanks = TRUE)
#' ex_data$outseq #view the output sequence, edges were reattached
#' #you will avoid data loss on edge of sequence, but errors in edge, or
#' #off target sequence will be present in the output
#' #
#' #build output sequence with only the COI-5P region
#' ex_data = outseq(ex_data, keep_flanks = FALSE)
#' ex_data$outseq #view the output sequence
#' #Ns added to the front to buffer trimmed region
#' #Note some sequence is lost due to the strange 
#' #path match that occurs at the front of the sequence.
#' ex_data$data$path
#' @name outseq
#' @export
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
  }else if(keep_flanks == 'right'){

    x$outseq = paste(c(x$adjusted_sequence, #the 'sweet spot' that is denoised
                         x$data$adjusted_trimmed, # part removed in adj_seq
                         x$frame_dat$removed_end, #part removed in second frame() call
                         x$data$raw_removed_end), # part removed from the end of first frame() call
                       collapse = "")
    
    if(!is.null(names(x$adjusted_sequence)) || !is.null(names(x$frame_dat$removed_lead))|| !is.null(names(x$data$raw_removed_front))){
      x$outphred = paste(c(names(x$adjusted_sequence), #the 'sweet spot' that is denoised
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