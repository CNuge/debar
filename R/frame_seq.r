#' build an DNAbin with ape.
#'
#' Switches the dashes in the seq - to n
#'
#' @keywords internal
individual_DNAbin = function(dna_string){
  return(ape::as.DNAbin(strsplit(gsub('-', 'n', as.character(tolower(dna_string))),"")))
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
  fwd_score = fwd_ntPHMMout[['score']]
  rev_score = rev_ntPHMMout[['score']]
  
  if(rev_score > fwd_score){
    outlist = list(raw = rev_x, ntBin = rev_ntBin, ntPHMMout = rev_ntPHMMout)
  }else{
    outlist = list(raw = x, ntBin = fwd_ntBin, ntPHMMout = fwd_ntPHMMout)    
  }
  return(outlist)
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
set_frame = function(org_seq_vec, path_out){
  
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
  
  if(!is.null(front) && !is.null(names(org_seq_vec))){
    names(front) = rep("~", length(front))
  }
  
  return(list(framed_seq = paste(c(front,org_seq_vec[org_seq_start:org_seq_end]),collapse= ""),
              framed_vec = c(front,org_seq_vec[org_seq_start:org_seq_end]),
              trimmed_seq = org_seq_vec[org_seq_start:org_seq_end],
              removed_lead = removed_lead,
              removed_end = removed_end,
              seq_start = org_seq_start,
              seq_end = org_seq_end,
              path_start = path_start,
              path_end = path_end,
              front = front))
}


#' Take a DNAseq object and isolate the COI-5P region. 
#' 
#' Raw DNA sequence is taken and compared against the nucleotide profile hidden Markov model (PHMM), 
#' generating the statistical information later used to apply corrections to the sequence read. The 
#' speed of this function and how it operates are dependent on the paramaters chosen. When dir_check == TRUE
#' the function will compare the forward and reverse compliment against the PHMM. Since the PHMM comparison is
#' a computationally intensive process, this option should be set to FALSE if you're data are all forward reads.
#' 
#' The min match and max inserts provide the minimal amount of matching to the PHMM requried for a sequence to not
#' be flagged for rejection. If you are dealing with sequence fragments much shorter than the entire COI-5P region you
#' should lower these values.
#' 
#' @param x a DNAseq class object.
#' @param dir_check A boolean indicating if both the forward and reverse compliments of a sequence should 
#' be checked against the PHMM. Default is TRUE.
#' @param double_pass A boolean indicating if a second pass through the Viterbi algorithm should be conducted for sequences
#' that had leading nucleotides not matching the PHMM. This improves the accurate establishment of reading frame and will
#' reduce false rejections by the amino acid check, but this comes at a cost of additional processing time. Default is TRUE.
#' @param min_match The minimum number of sequential matches to the PHMM for a sequence to be denoised.
#' Otherwise flag the sequence as a reject.
#' @param max_inserts The maximum number of sequention insert states occuring in a sequence 
#' (including the flanking regions). If this number is exceeded than the entire read will be labelled for rejection. 
#' Default is 400.
#' @param terminate_rejects Should a check be made to enusre minimum homology of the input sequence to the PHMM.
#' Makes sure the match conditions are met prior to continuing with the framing of the sequence. If conditions not
#' met then the function is stopped and the sequence labelled for rejection.
#' @param ... additional arguments to be passed between methods.
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @examples
#' #previously called
#' hi_phred = paste0(rep("~~~", 217), collapse="")
#' dat1 = DNAseq(example_nt_string_errors, name = "err_seq1", phred = hi_phred)
#' dat1= frame(dat1)
#' dat1$frame_dat # a labelled list with framing information is produced
#' dat1$reject == FALSE # the match criteria were met, not labelled for rejection
#' dat1$data$path #one can call and view the path diagram produced by the PHMM comparison.
#' @export
#' @name frame
frame = function(x, ...){
  UseMethod("frame")
} 


#' @rdname frame
#' @export
frame.DNAseq = function(x, ..., dir_check = TRUE, double_pass = TRUE, min_match = 100, max_inserts = 400, terminate_rejects = TRUE){
  
  if(dir_check == TRUE){
    dir_out = dir_check(x$raw)
    #parse the outputs from the optimal direction.
    x$raw = dir_out$raw
    x$data$ntBin = dir_out$ntBin
    x$data$ntPHMMout = dir_out$ntPHMMout
  }else{
    x$data$ntBin = individual_DNAbin(toupper(x$raw))
    x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
  }
  
  #take the optimal direction's output string and reject if:
  #a. There is a large string of 2s (indicative of a chimera) 
  #OR
  #b. There is not a continious match to the PHMM of at least the designated min_match length.
  if(grepl(paste(rep("2", max_inserts), collapse = ""), paste( x$data$ntPHMMout[['path']], collapse = "")) ||
     !grepl(paste(rep("1", min_match), collapse = ""), paste( x$data$ntPHMMout[['path']], collapse = "")) ){
    
    x$reject = TRUE
  }else{
    x$reject = FALSE
  }
  #if the reject condition is met and we are early terminating rejects,
  #then clean the object and return as opposed to establishing reading frame.
  if(terminate_rejects == TRUE && x$reject == TRUE){
    x$data$path = x$data$ntPHMMout[['path']]
    
    x$data$ntBin = NULL
    x$data$ntPHMMout = NULL
    
    return(x) 
  }
    
  #turn the raw string into a vector, add phred labels as well.
  org_seq_vec = strsplit(tolower(x$raw), split='')[[1]]
  #add the labels
  if(!is.null(x$phred)){
    names(org_seq_vec) = strsplit(x$phred, split='')[[1]]
  }
  
  #check for leading inserts, if present remove them and reframe the sequence for higher accuracy.
  if(leading_ins(x$data$ntPHMMout[['path']]) && double_pass == TRUE){
    temp_frame = set_frame(org_seq_vec, x$data$ntPHMMout[['path']]) #run initial framing of the full sequence
    x$data$raw_removed_front = temp_frame[['removed_lead']] # keep any tirmmed edges of the raw sequence separate
    x$data$raw_removed_end = temp_frame[['removed_end']]    #
    x$data$len_first_front = length(temp_frame$front)
    trim_temp = temp_frame[['framed_seq']]
    x$data$ntBin = individual_DNAbin(trim_temp) #rerun the trimmed sequence through the PHMM, to more accurately establish frame
    x$data$ntPHMMout = aphid::Viterbi(nt_PHMM, x$data$ntBin, odds = FALSE)
    
    #set the reading frame for the final sequence, produces data to
    #constrain the operation of the adjust seqence function
    x$frame_dat = set_frame(temp_frame$framed_vec, x$data$ntPHMMout[['path']])
    x$data$path = x$data$ntPHMMout[['path']]
  }else{
    x$data$raw_removed_front = c()
    x$data$raw_removed_end = c()
    #set the reading frame for the final sequence, produces data to
    #constrain the operation of the adjust seqence function
    x$frame_dat = set_frame(org_seq_vec, x$data$ntPHMMout[['path']])
    x$data$path = x$data$ntPHMMout[['path']]
  }

  #remove the aphid and ape strucutres to minimize the memory footprint.
  x$data$ntBin = NULL
  x$data$ntPHMMout = NULL
  
  return(x)
}
