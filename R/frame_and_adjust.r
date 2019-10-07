#' build an DNAbin with ape.
#'
#' Switches the dashes in the seq - to n
#'
#' @keywords internal
individual_DNAbin = function(dna_string){
  return(ape::as.DNAbin(strsplit(gsub('-', 'n', as.character(tolower(dna_string))),"")))
}

#' Check for a large number of leading inserted bases,
#' if this is the case, TRUE is returned and the PHMM
#' should be run a second time on the truncated data.
#' @keywords internal
leading_ins = function(seq_path){
  #iterate along the sequence, if we hit 5 insertions
  #before 5 matches then there is a problem
  matches = 0
  ins = 0
  for(x in seq_path){
    if (x == 1){
      matches = matches + 1
    } else if (x == 2){
      ins = ins + 1
    }
    
    if (matches == 5){
      return(FALSE)
    }
    if (ins == 5){
      return(TRUE)
    }
  }
}

#' check sequence for an early large string of deletions, if it exists then
#' return the starting index by which to slice the path and the string
#' @keywords internal
ins_front_trim = function(path_out, search_scope = 15){
  if(sum(path_out[1:search_scope] == 2)>2){
    run_pos = min(which(path_out == 2))
    for(i in run_pos:length(path_out)){
      if(path_out[i] != 2){
        return(i)
      }
    }
  }
  return(FALSE)
}


#' Take an input sequence and get it into the reading frame.
#' Uses the path of the sequence through the PHMM to located
#' the first long contigious stretch of 1s (matches)
#' Note - there is an imbalance in operation, namely the function
#' checks for a run of 10 additional 1s on the front of the sequence, and for
#' a run of only 3 additional matches on the back. The front of the sequence is high
#' priority because it is necessary to ensure the sequence is in frame.
#' whereas matching on the back is more lenient as insertions and deletions
#' can be tolerated here without large implications for the rest of the sequence
#' @keywords internal
set_frame = function(org_seq , path_out){
  org_seq_vec = strsplit(tolower(org_seq), split='')[[1]]
  
  adj_for_dels = ins_front_trim(path_out)
  front = c()
  removed_lead = c()
  
  #If there are a large number (>2) of deletes in the first 15bp
  #then adjust the sequence so the end of the delete run is the starting point
  
  
  org_seq_start = 1
  org_seq_end = length(org_seq_vec)
  
  path_start = 1
  path_end = length(path_out)
  
  if(adj_for_dels != FALSE){
    org_seq_start = adj_for_dels
    path_start = adj_for_dels
  }
  
  
  for( i in 1:length(path_out) ){
    #0 = D
    if( path_out[i] == 0 ){
      #there was a bp missing in the original seq
      #add a placeholder to the new seq, don't advance on original pointer
      front = c(front, '-')
      #2 = I
    }else if( path_out[i] == 2 ){
      #there was an extra bp at the front of the sequence
      #not represented in the PHMMs, skip this bp in the original seq
      org_seq_start = org_seq_start + 1
      #for the first match seen, check to make sure it isn't a single match 
      #or codon of dangling matches in a sea of inserts or deletes
    }else if( path_out[i] == 1 ){
      if (2 %in% path_out[(i+1):(i+10)] ){
        org_seq_start = org_seq_start + 1
      }else if (0 %in% path_out[(i+1):(i+10)] ){
        org_seq_start = org_seq_start + 1
        front = c(front, '-')
      }else{
        path_start = i
        break
      }
    }
  }
  
  #trim the back of the sequence
  for( i in length(path_out):1 ){
    #0 = D
    if( path_out[i] == 0 ){
      #there was a bp missing at the back of the original seq
      #just continue to the next position
      next
      #2 = I
    }else if( path_out[i] == 2 ){
      #there is an extra base pair at the end of the original sequence
      #not represented in the PHMMs, remove this trailing end of the original sequence
      org_seq_end = org_seq_end - 1
    }else if ( path_out[i] == 1){
      #in either of these instances we want to trim the last bp, as its
      #dangling in a sea of inserts or deletes and likely a random profile match
      #if (path_out[i-1] == 2 & path_out[i-2] == 2){
      if (2 %in% path_out[(i-3):(i-1)]){
        org_seq_end = org_seq_end - 1
      }else if (0 %in% path_out[(i-3):(i-1)]){
        org_seq_end = org_seq_end - 1
      }else{
        path_end = i
        break
      }
    }
  }
  

  if(org_seq_start != 1){
    removed_lead = org_seq_vec[1:(org_seq_start-1)]
  }else{
    removed_lead = c()
  }
  
  if(org_seq_end < length(org_seq_vec)){
    removed_end = org_seq_vec[(org_seq_end+1):length(org_seq_vec)]
  }else{
    removed_end = c()
  }
  
  return(list(framed_seq = paste(c(front,org_seq_vec[org_seq_start:org_seq_end]),collapse= ""),
              trimmed_seq = org_seq_vec[org_seq_start:org_seq_end],
              removed_lead = removed_lead,
              removed_end = removed_end,
              seq_start = org_seq_start,
              seq_end = org_seq_end,
              path_start = path_start,
              path_end = path_end,
              front = front))
}

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

#' Adjust the DNA sequence based on the ntPHMM path
#'
#' @param frame_dat The DNAseq's framing data - generated by set_frame()
#' @param path_out The nucleotide PHMM path for the sequence.
#' @param censor_length number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters 
#' @keywords internal
adj_seq = function(frame_dat, path_out, censor_length = 3){
  
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

  #adjustments to skip because they occur as a codon
  triple_inserts = triple_ins(path_out, path_pos, path_end)

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
          new_seq = c(new_seq, '-')
          
          new_pos = new_pos + 1
          
          censor_0s = c(censor_0s, new_pos)
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
        #there was a bp insertion, skip this bp in the original seq
        org_seq_pos = org_seq_pos + 1
        
        censor_2s = c(censor_2s, new_pos)
      }
    }

  }
  
  #if the end of the original seq end was not reached, add the remained to the 
  
  if(censor_length != 0){

    #censor_0s
    if(!is.null(censor_0s)){
      censor_0s_masks = unlist(lapply(censor_0s, function(pos){
        (pos-censor_length):(pos+censor_length)
      }))		
      censor_0s_masks = unique(censor_0s_masks)
      censor_0s_masks = censor_0s_masks[censor_0s_masks>0]
      censor_0s_masks = censor_0s_masks[censor_0s_masks<=length(new_seq)]
      
      new_seq[censor_0s_masks] = "-"
      
    }
    #censor_2s
    if(!is.null(censor_2s)){
      censor_2s_masks = unlist(lapply(censor_2s, function(pos){
        (pos-(censor_length-1)):(pos+censor_length)
      }))
      censor_2s_masks = unique(censor_2s_masks)
      censor_2s_masks = censor_2s_masks[censor_2s_masks>0]
      censor_2s_masks = censor_2s_masks[censor_2s_masks<=length(new_seq)]
      
      new_seq[censor_2s_masks] = "-"		
    }
  }

  adj_count = length(c(censor_0s,censor_2s))
  
  return(list(c(frame_dat$front, new_seq), adj_count))
}

#' Take a DNAseq object and put sequences in common reading frame
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @examples
#' #previously called
#' ex_data = DNAseq(ex_ccs_read_list, name = 'ex1')
#' #frame the ccs reads, will utilize the PHMM corresponding to the order specified in the object.
#' ex_data =  frame(ex_data)
#' @export
#' @name frame
frame = function(x, ...){
  UseMethod("frame")
} 


#' @rdname frame
#' @export
frame.DNAseq = function(x, ...){
  
  x$data$ntBin = individual_DNAbin(toupper(x$raw))
  
  x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
  
  if(leading_ins(x$data$ntPHMMout[['path']])){
    temp_frame = set_frame(x$raw, x$data$ntPHMMout[['path']])
    x$data$raw_removed_front = temp_frame[['removed_lead']]
    x$data$raw_removed_end = temp_frame[['removed_end']]
    x$data$len_first_front = length(temp_frame$front)
    trim_temp = temp_frame[['framed_seq']]
    x$data$ntBin = individual_DNAbin(trim_temp)
    x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
  }else{
    x$data$raw_removed_front = c()
    x$data$raw_removed_end = c()
    trim_temp = x$raw
  }
  x$frame_dat = set_frame(trim_temp, x$data$ntPHMMout[['path']])
  x$data$path = x$data$ntPHMMout[['path']]
 
  x$data$ntBin = NULL
  x$data$ntPHMMout = NULL
  
  return(x)
}


#' Adjust the sequences based on the nt path outputs.
#'
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param censor_length the number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 3.
#' @return a class object of code{"ccs_reads"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @examples
#' #' #previously called
#' ex_data = DNAseq(example_nt_string, name = 'SSGBC787-14')
#' ex_data =  frame(ex_data)
#' #adjust the sequence with default censor length is 3
#' ex_data = adjust(ex_data)
#' #with a custom censorship size
#' ex_data = adjust(ex_data, censor_length = 5)
#' @export
#' @name adjust
adjust = function(x, ...){
  UseMethod("adjust")
} 

#' @rdname adjust
#' @export
adjust.DNAseq = function(x, ..., censor_length = 3){
  
  adj_out = adj_seq(x$frame_dat, x$data$path, censor_length = censor_length)

  x$adjusted_sequence = adj_out[[1]]
  x$adjustment_count = adj_out[[2]]
  return(x)
}
