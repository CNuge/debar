
#' Return a labelled list 
#' @keywords internal
phred2numeric = function(phred_string){
  outvec = gtools::asc(unlist(strsplit(phred_string, "")))
  return(outvec - 33)
}


#' Check the numeric phred scores and 
#'
#' @param x a DNAseq class object.
#' @param min_avg_qv
#' @param max_perc_low The maximum percentage of nucleotides in the string with QV values lower than 20. Default is 25%
#' @param max_perc_ultra_low The maximum percentage of nucleotides in the string with QV values lower than 10. Default is 5%
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
