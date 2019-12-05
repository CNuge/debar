#' Print a concise summary of the DNAseq object
#' @keywords internal
print.DNAseq = function(x, ...){
  desc_string = "A DNAseq object.\n"
  
  lines = c(desc_string)
  if(length(x$name) != 0){
    l2 = paste("Sample ID: ", x$name, "\n", sep = "")  
    lines = c(lines, l2)
  }
  l3 = paste("Raw Sequence:\n", 
             substr(x$raw,1, 25), "..." , substr(x$raw, (nchar(x$raw)-24), nchar(x$raw)), "\n",
             sep ="")
  lines = c(lines, l3)
  
  if(!is.null(x$adjustment_count)){
    l4 = paste("Corrections applied: ", x$adjustment_count , "\n",sep = "")
    lines = c(lines, l4)
    
  }
  
  if(!is.null(x$outseq)){
    l5 = paste("Output Sequence:\n",
               substr(x$outseq,1, 25), "..." , substr(x$outseq, (nchar(x$outseq)-24), nchar(x$outseq)), "\n",
               sep ="")
    lines = c(lines, l5)
  }
  
  cat(lines, sep="")
} 
