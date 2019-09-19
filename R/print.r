#' Print a concise summary of the ccs_reads object
#' @keywords internal
print.ccs_reads = function(x, ...){
  desc_string = paste(
    "A ccs_reads object comprised of ", length(x$sequences), " unique ccs reads.\n",
    "Taxonomic order: ", x$order, 
    sep="")
  
  lines = c(desc_string)
  if(length(x$id) != 0){
    l2 = paste("\nSample ID: ", x$id, '\n', sep = "")  
    lines = c(lines, l2)
  }
  cat(lines, sep="")
}
# 
