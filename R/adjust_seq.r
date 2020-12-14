
#' Look for triple inserts in the PHMM path.
#' @keywords internal
triple_ins = function(x, path_start, path_end){
  inserts = c(path_start:path_end)[x[path_start:path_end] == 0]
  triples = c()
  if(length(inserts) > 2){
    for(i in 1:(length(inserts)-2)){
      potential = inserts[i]:(inserts[i]+2) 
      observed = inserts[i:(i+2)]
      if(isTRUE(all.equal(potential, observed))){
        triples = c(triples, observed)
      }
    }
    return(triples)
  }else{
    return(c())
  }
}

#' Look for triple deletes in the PHMM path.
#' @keywords internal
triple_dels = function(x, path_start, path_end){
  deletes = c(path_start:path_end)[x[path_start:path_end] == 2]
  triples = c()
  if(length(deletes) > 2){
    for(i in 1:(length(deletes)-2)){
      potential = deletes[i]:(deletes[i]+2) 
      observed = deletes[i:(i+2)]
      if(isTRUE(all.equal(potential, observed))){
        triples = c(triples, observed)
      }
    }
    return(triples)
  }else{
    return(c())
  }
}


#' Adjust the DNA sequence based on the ntPHMM path.
#'
#' @param frame_dat The DNAseq's framing data - generated by set_frame()
#' @param path_out The nucleotide PHMM path for the sequence.
#' @param censor_length Number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters.
#' @param added_phred The phred character to use for characters inserted into the original sequence.
#' @return A named list.
#' @keywords internal
adj_seq = function(frame_dat, path_out, censor_length = 3, added_phred = "*"){
  
  org_seq_vec = frame_dat$trimmed_seq
  #or pass in the bin value and use:
  #org_seq_vec = as.character(org_seq)[[1]]
  
  # build a new sequence from scratch
  #note: the front of the sequence (leading dashes)
  #is not added until the last line, after adjustment and censorship
  #this is to prevent the dashes from shifting the frame and causing the 
  #positions of all of the corrections to have to be adjusted
  new_seq = c()
  #start at the first positon in the path
  org_seq_pos = 1
  path_pos = frame_dat$path_start
  path_end = frame_dat$path_end
  new_pos = 1
  
  censor_0s = c() #these were deleted from true seq  
  censor_2s = c() #these were inserted into true seq
  match_seen = FALSE
  first_0 = TRUE

  #adjustments to skip because they occur as a codon
  triple_inserts = triple_ins(path_out, path_pos, path_end)
  triple_deletes = triple_dels(path_out, path_pos, path_end)
  
  for(i in path_pos:path_end){
    #this prevents additions beyond the end of the original sequence.
    if(org_seq_pos > length(org_seq_vec)){
      break
    }
    #0 = D
    if(path_out[i] == 0){
      #there was a bp missing in the original seq
      #add a placeholder to the new seq, don't advance on original pointer
      
      #NOTE: this is to truncate the sequence if a large string of
      #gap characters is encountered, this prevents incorrect additions
      #The information is not discarded, but rather added to an
      #additional part of the output vector, so that it can be recovered if
      #the edges of the sequences are not being removed
      if(paste(path_out[i:(i+4)], collapse='') == "00000"){
        break        
      }
      #NOTE: this if clause is to prevent addition of an extra placeholder
      #in cases where adjust seq has already added a placeholder
      if(match_seen == TRUE){
        #only make the placeholder addition if the 
        #the 0 in the path is not found as a member
        #of a triple, these are missing codons and likely true
        if(!i %in% triple_inserts){
          
          if(first_0 == TRUE){
            #First add the bp at this position normally if previous position not a 0
            new_seq = c(new_seq, org_seq_vec[org_seq_pos])
            org_seq_pos = org_seq_pos + 1
            
            new_pos = new_pos + 1
            match_seen = TRUE
            first_0 = FALSE
          }
          #then add a placeholder behind it
          add_char = '-'
          
          if(!is.null(names(org_seq_vec))){
            names(add_char) = added_phred
          }
          
          new_seq = c(new_seq, add_char)
          new_pos = new_pos + 1
          censor_0s = c(censor_0s, new_pos)
          if(i < path_end){
            if(path_out[(i+1)] != 0){
              first_0 = TRUE
            }
          }
        }
      }
      
      #1 = M
    }else if(path_out[i] == 1){
      new_seq = c(new_seq, org_seq_vec[org_seq_pos])
      org_seq_pos = org_seq_pos + 1
      
      new_pos = new_pos + 1
      match_seen = TRUE
      
      #2 = I
    }else if(path_out[i] == 2){
      if(match_seen == TRUE){
        #if there is a large run of 2s then cease the adjustment loop, sequence if
        #off the rails and should not be adjusted further
        if(paste(path_out[i:(i+4)], collapse='') == "22222"){
          break        
        }
        if(!i %in% triple_deletes){
          #there was a bp insertion, skip this bp in the original seq
          org_seq_pos = org_seq_pos + 1
          
          censor_2s = c(censor_2s, new_pos)
        }
      }
    }
  }
  
  if(org_seq_pos < length(org_seq_vec)){
    trimmed_end = org_seq_vec[org_seq_pos:length(org_seq_vec)]
  }else{
    trimmed_end = c()
  }
    
  if(censor_length != 0){

    #censor_0s
    if(!is.null(censor_0s)){
      censor_0s_masks = unlist(lapply(censor_0s, function(pos){
        (pos-censor_length-1):(pos+censor_length-1)
      }))		
      censor_0s_masks = unique(censor_0s_masks)
      censor_0s_masks = censor_0s_masks[censor_0s_masks>0]
      censor_0s_masks = censor_0s_masks[censor_0s_masks<=length(new_seq)]
      
      new_seq[censor_0s_masks] = "-"
      
    }
    #censor_2s
    if(!is.null(censor_2s)){
      censor_2s_masks = unlist(lapply(censor_2s, function(pos){
        (pos-(censor_length)):(pos+censor_length-1)
      }))
      censor_2s_masks = unique(censor_2s_masks)
      censor_2s_masks = censor_2s_masks[censor_2s_masks>0]
      censor_2s_masks = censor_2s_masks[censor_2s_masks<=length(new_seq)]
      
      new_seq[censor_2s_masks] = "-"		
    }
  }

  adj_count = length(c(censor_0s,censor_2s))
  
  return(list(c(frame_dat$front, new_seq), adj_count, trimmed_end, censor_0s, censor_2s))
}

#' Adjust the sequences based on the nt path outputs.
#' 
#' Based on the PHMM path generated by the frame function, the sequence is the adjusted. Adjustments are limited
#' to the 657bp region represeneted by the PHMM (and the part of the input sequence matching this region). Censorship
#' can be applied around the corrections. This limits the number of indel errors missed by the PHMM correction algorithm,
#' but comes at a cost of lost DNA sequence. The default of 7 is a conservative paramater meant to lead to the coverage of
#' greater than 95% of indel errors.
#' 
#' If the DNAseq object contains PHRED scores, the PHRED string will be adjusted along with the DNA sequence (corresponding)
#' value removed when a bp removed. The 'added_phred' value indicated the phred chracter to be added to the string when a 
#' placeholder nucleotide is added to the string to account for a deletion. Default is "*" which indicates a score of 9 (a 
#' relitavely low quality base).
#'
#' @param x A DNAseq class object.
#' @param ... Additional arguments to be passed between methods.
#' @param censor_length The number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 7.
#' @param added_phred The phred character to use for characters inserted into the original sequence. Default is "*".
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @examples
#' #previously called
#' ex_data = DNAseq(example_nt_string_errors, name = 'error_adj_example')
#' ex_data =  frame(ex_data)
#' #adjust the sequence with default censor length is 7
#' ex_data = adjust(ex_data)
#' ex_data$adjusted_sequence #output is a vector, use outseq to build the string
#' #with a custom censorship size
#' ex_data = adjust(ex_data, censor_length = 5)
#' ex_data$adjusted_sequence #less flanking base pairs turned to placeholders
#' ex_data$adjustment_count #get a count of the number of adjustments applied
#' @export
#' @name adjust
adjust = function(x, ...){
  UseMethod("adjust")
} 

#' @rdname adjust
#' @export
adjust.DNAseq = function(x, ..., censor_length = 7,  added_phred = "*"){
  
  adj_out = adj_seq(x$frame_dat, x$data$path, censor_length = censor_length, added_phred = added_phred)

  x$adjusted_sequence = adj_out[[1]]
  x$adjustment_count = adj_out[[2]]
  x$data$adjusted_trimmed = adj_out[[3]]
  x$data$deletes_changed = adj_out[[4]]
  x$data$inserts_changed = adj_out[[5]]
  return(x)
}

