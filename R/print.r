#' Print a concise summary of the DNAseq object
#' @keywords internal
print.DNAseq = function(x, ...){
  desc_string = "A DNAseq object.\n"
  
  lines = c(desc_string)
  if(length(x$id) != 0){
    l2 = paste("\nSample ID: ", x$id, "\n", sep = "")  
    lines = c(lines, l2)
  }
  l3 = paste("Raw Sequence\n", 
             substr(x$raw,1, 25), "..." , substr(x$raw, (nchar(x$raw)-24), nchar(x$raw)), "\n",
             sep ="")
  lines = c(lines, l3)
  
  #if()
  
  cat(lines, sep="")
} 
