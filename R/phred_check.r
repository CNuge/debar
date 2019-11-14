
#' Return a labelled list 
#' @keywords internal
phred2numeric = function(phred_string){
  outvec = gtools::asc(unlist(strsplit(phred_string, "")))
  return(outvec - 33)
}


#' Check that numeric phred scores indicate the sequence meets quality standards.
#'
#' This function is used to separate out low quality reads from a data set. Sequences with phred values not meeting
#' the minimum quality standards specified. In the background it is converting the string of phred scores to numeric
#' values and then comparing these against the thresholds to make sure the quality standards are met. If not, the
#' returned DNAseq class will have a reject field equal to true.
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param min_avg_qv The minimum average phred score for a read to be retained.
#' @param max_perc_low The maximum frequency of nucleotides in the string with QV values lower than 20. Default is 0.25
#' @param max_perc_ultra_low The maximum frequency of nucleotides in the string with QV values lower than 10. Default is 5
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @examples
#' # simulate a high quality phred score
#' hi_phred = paste0(rep("~~~", 217), collapse="")
#' dat1 = DNAseq(example_nt_string_errors, name = "phred_good", phred = hi_phred)
#' dat1 = phred_check(dat1) #default params
#' #check if the read should be rejected
#' dat1$reject == FALSE
#' # simulate a low quality phred score
#' low_phred = paste0(rep("*0-", 217), collapse="")
#' dat2 = DNAseq(example_nt_string_errors, name = "phred_bad", phred = low_phred)
#' dat2 = phred_check(dat2) #default params
#' #check if the read should be rejected
#' dat2$reject == TRUE
#' @export
#' @name phred_check
phred_check = function(x, ...){
  UseMethod("phred_check")
}

#' @rdname phred_check
#' @export
phred_check.DNAseq = function(x, ..., 
                              min_avg_qv = 20,
                              max_perc_low = 0.25,
                              max_perc_ultra_low = 0.05){

  if(is.null(x$phred)){
    return(x)
  }
  num_phred = phred2numeric(x$phred)
  names(num_phred) = NULL
  avg_phred = mean(num_phred)
  
  if(avg_phred < min_avg_qv){
    x$reject = TRUE
    return(x)
  }
  
  num_low = length(num_phred[num_phred < 20])
  if(num_low > (length(num_phred) * max_perc_low)){
    x$reject = TRUE
    return(x)
  }
  
  num_ultra_low = length(num_phred[num_phred < 10])
  if(num_ultra_low > (length(num_phred) * max_perc_ultra_low)){
    x$reject = TRUE
    return(x)
  }
  
  x$reject = FALSE
  x
}
