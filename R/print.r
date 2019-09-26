#' Print a concise summary of the DNAseq object
#' @keywords internal
print.DNAseq = function(x, ...){
  desc_string = "A DNAseq object.\n"
  
  lines = c(desc_string)
  if(length(x$id) != 0){
    l2 = paste("\nSample ID: ", x$id, "\n", sep = "")  
    lines = c(lines, l2)
  }
  l3 = paste("Sequence\n", x$sequence, "\n",
             sep ="")
  lines = c(lines, l3)
  cat(lines, sep="")
} 
