
#' Read in raw data from a fastq file.
#' 
#' 
#' @param x The name of the fastq file to read data from.
#' @param keep_quality Boolean indicating if the Phred quality scores should be 
#' retained in the output dataframe. Default is TRUE
#' @export
#' @name read_fastq
read_fastq = function(x, keep_quality = TRUE){
  
  records = list()  
  
  n = 4
  lines = c()
  
  for(i in readLines(x)){
    
    lines = c(lines, i)
    
    if(length(lines) == n){
      records[[length(records)+1]] =  c(substr(lines[1], 2, nchar(lines[1])),
                                      lines[2],
                                      lines[4])
      lines = c()
    }
  }
  
  out = data.frame(do.call(rbind, records) , stringsAsFactors = FALSE)
  names(out)= c("header_data", "sequence", "quality")
  
  if(keep_quality == FALSE){
    out$quality = NULL
  }
  
  return(out)
}


#' Read in raw data from a fasta file.
#' 
#' @param x The name of the fasta file to read data from
#' @export
#' @name read_fasta
read_fasta = function(x){
  
  data = readLines(x)
  
  head_line = data[seq(1,length(data), 2)]
  head_line = substr(head_line, 2, nchar(head_line))
  
  records = data.frame(header = head_line,
                       sequence = data[seq(2,length(data), 2)], 
                       stringsAsFactors = FALSE)
  
  return(records)
}


#' Output the denoised consensus sequence to a fasta file.
#'
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param filename The name of the file to output the data to. Default is "denoised.fasta".
#' @param append Should the ccs consensus sequence be appended to the output file?(TRUE) 
#' Or overwrite the file?(FALSE) Default is TRUE.
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @export
#' @name write_fasta
write_fasta = function(x, ...){
  UseMethod("write_fasta")
  
}

#' @rdname write_fasta
#' @export
write_fasta.DNAseq = function(x, ...,
                                 filename = "denoised.fasta", 
                                 append = TRUE){

  outstring = paste(">", x$name, "\n",
                    x$outseq, sep = '')
  
  write(outstring, file = filename, append = append)
}


#' Output the denoised sequence to a fastq format with placeholder phred scores.
#'
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param filename The name of the file to output the data to. Default is "denoised.fasta".
#' @param append Should the ccs consensus sequence be appended to the output file?(TRUE) 
#' Or overwrite the file?(FALSE) Default is TRUE.
#' @param keep_phred Should the original PHRED scores be kept in the output? Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'.
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @export
#' @name write_fastq
write_fastq = function(x, ...){
  UseMethod("write_fastq")
}


#' @rdname write_fastq
#' @export
write_fastq.DNAseq = function(x, ...,
                                 filename = "denoised.fastq", 
                                 append = TRUE, 
                                 keep_phred = TRUE,
                                 phred_placeholder = "#"){

  if(keep_phred == TRUE){
    if(is.null(x$phred)){
      stop("Cannot keep the phred scores for a DNAseq with no phred inputs")
    }
    outstring = paste("@", x$name, "\n",
                      x$outseq, "\n",
                      "+\n",
                      x$outphred, sep="")
  }else{
    outstring = paste("@", x$name, "\n",
                      x$outseq, "\n",
                      "+\n",
                      paste(rep(phred_placeholder, times = nchar(x$outseq)), collapse = ""), sep="")
  }
  
  write(outstring, file =  filename, append = append)
}

#' A wrapper function to deploy the fastq and fata output functions.
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param filename The name of the file to output the data to. Default is "denoised.fasta".
#' @param outformat The format of the output data, either fasta for fastq. Default is fastq.
#' @param append Should the ccs consensus sequence be appended to the output file?(TRUE) 
#' Or overwrite the file?(FALSE) Default is TRUE.
#' @param keep_phred Should the original PHRED scores be kept in the output? Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'.
#' @return a class object of code{"DNAseq"}
#'
#'
write_wrapper = function(x, ...){
  UseMethod("write_wrapper")
}

#' @rdname write_wrapper
#' @export
write_wrapper.DNAseq = function(x, ...,
                              filename = "denoised.fastq", 
                              outformat = "fastq",
                              append = TRUE, 
                              keep_phred = TRUE,
                              phred_placeholder = "#"){
  if(outformat != "fastq" && outformat != "fasta" && outformat != "none"){
    stop("Invalid output format! Must be one of: 'fasta', 'fastq' or 'none'")
  }
  
  if(outformat == "fasta"){
    if(is.null(filename)){
      write_fasta(x, append = append, ...)
    }else{
      write_fasta(x, filename = filename, append = append, ...)
    }
  }
  
  if(outformat == "fastq"){
    if(is.null(filename)){
      write_fastq(x, keep_phred = keep_phred, 
                  phred_placeholder = phred_placeholder, 
                  append = append, ...)
    }else{
      write_fastq(x, keep_phred = keep_phred,
                  filename = filename, 
                  phred_placeholder = phred_placeholder, 
                  append = append, ...)
    }
  }
  x
}