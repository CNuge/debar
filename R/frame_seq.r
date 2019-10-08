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
  
  
  for( i in path_start:length(path_out) ){
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
      if (2 %in% path_out[(i-10):(i-1)]){
        org_seq_end = org_seq_end - 1
      }else if (0 %in% path_out[(i-10):(i-1)]){
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




#' Return the reverse compliment for a DNA sequence
#' @keywords internal
rev_comp = function(x){
  return(seqinr::c2s(rev(seqinr::comp(seqinr::s2c(x)))))
}


#' Take an input sequence and align both the forward and reverse compliments to the PHMM
#' 
#' The fuction returns the sequence, DNAbin and Path and score of the optimal orientation.
#' Optimal orientation is determined by the direction with the longer string of consecutive 
#' ones in the path
#' @param x a DNAseq class object.
#' 
dir_check = function(x){
  #run this for fwd
  fwd_ntBin = individual_DNAbin(toupper(x))
  fwd_ntPHMMout = aphid::Viterbi(nt_PHMM, fwd_ntBin, odds = FALSE)
  
  rev_x = rev_comp(x)
  rev_ntBin = individual_DNAbin(toupper(rev_x))
  rev_ntPHMMout = aphid::Viterbi(nt_PHMM, rev_ntBin, odds = FALSE)
  
  #compare the PHMM scores and return the one with the higher logL
  
  
}

#' Take a DNAseq object and put sequences in common reading frame
#' @param x a DNAseq class object.
#' @param dir_check Should both the forward and reverse compliments be considered?
#' @param min_logL The minimum log likelihood value for a sequence to be denoised.
#' 
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
  
  if(dir_check == TRUE){
    dir_out = dir_check(x$raw)
    #parse the outputs from the optimal direction.
    x$data$ntBin = 
    x$data$ntPHMMout = 
  }else{
    x$data$ntBin = individual_DNAbin(toupper(x$raw))
    x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
  }
  
  #TODO - finish this check of standards.
  #if the score of the PHMMout
  if(x$data$ntPHMMout < min_logL){
    #add a signal for the denoise function to skip additional steps
    #break and don't waste time with the rest.

  }
  
  #check for leading inserts, if present remove them and reframe the sequence for higher accuracy.
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
  #set the reading frame for the final sequence, produces data to
  #constrain the operation of the adjust seqence function
  x$frame_dat = set_frame(trim_temp, x$data$ntPHMMout[['path']])
  x$data$path = x$data$ntPHMMout[['path']]
 
  #remove the aphid and ape strucutres to minimuze memory footprint.
  x$data$ntBin = NULL
  x$data$ntPHMMout = NULL
  
  return(x)
}
