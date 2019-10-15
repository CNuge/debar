
#' Read in raw data from a fastq file.
#' 
#' 
#' @param x The name of the fastq file to read data from.
#' @param keep_quality Boolean indicating if the Phred quality scores should be 
#' retained in the output dataframe. Default is TRUE
#' @examples
#' filename = system.file('extdata/ccs_subset.fastq', package = 'seqdenoise')
#' data = read_fastq(filename)
#' @export
#' @name read_fastq
read_fastq = function(x, keep_quality = TRUE){
  
  if(keep_quality == TRUE){
    records = data.frame(header_data = character(),
                         sequence = character(), 
                         quality = character(),
                         stringsAsFactors = FALSE)
  }else{
    records = data.frame(header_data = character(),
                         sequence = character(), 
                         stringsAsFactors = FALSE)
  }
  
  n = 4
  lines = c()
  
  for(i in readLines(x)){
  
    lines = c(lines, i)
    
    if(length(lines) == n){
      if(keep_quality == TRUE){
        records = rbind(records, data.frame(header_data = lines[1], 
                                            sequence = lines[2], 
                                            quality = lines[4], 
                                            stringsAsFactors = FALSE))
      }else{
        records = rbind(records, data.frame(header_data = lines[1], 
                                            sequence = lines[2], 
                                            stringsAsFactors = FALSE))
      }
      lines = c()
    }  
  }
  
  return(records)
}



#' Read in raw data from a fasta file.
#' 
#' @param x The name of the fasta file to read data from
#' @examples
#' filename = system.file('extdata/ccs_subset.fasta', package = 'coiDenoiser')
#' data = read_fasta(filename)
#' @export
#' @name read_fasta
read_fasta = function(x){
  
  data = readLines(x)
  
  records = data.frame(header = data[seq(1,length(data), 2)],
                       sequence = data[seq(2,length(data), 2)], 
                       stringsAsFactors = FALSE)
  
  return(records)
}


#' Output the denoised consensus sequence to a fasta file.
#'
#' 
#' @param x a DNAseq class object.
#' @param ... additional arguments to be passed between methods.
#' @param keep_flanks Default is TRUE.
#' @param filename The name of the file to output the data to. Default is "denoised.fasta".
#' @param append Should the ccs consensus sequence be appended to the output file?(TRUE) 
#' Or overwrite the file?(FALSE) Default is TRUE.
#' @return a class object of code{"DNAseq"} 
#' @seealso \code{\link{DNAseq}}
#' @seealso \code{\link{frame}}
#' @seealso \code{\link{adjust}}
#' @examples
#' #previously called
#' ex_data = DNAseq(example_nt_string, name = 'ex1')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data)
#' ex_data = consensus(ex_data)
#' #write to a fasta file with the default output file
#' write_fasta(ex_data)
#' #specify the path to a custom output file, overwrite its current contents
#' write_fasta(ex_data, filename = 'example_out.fasta', append = FALSE)
#' @export
#' @name write_fasta
write_fasta = function(x, ...){
  UseMethod("write_fasta")
  
}

#' @rdname write_fasta
#' @export
write_fasta.DNAseq = function(x, ...,
                                 keep_flanks = TRUE,       
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
#' @examples
#' #previously called
#' ex_data = DNAseq(example_nt_string, name = 'ex1')
#' ex_data =  frame(ex_data)
#' ex_data = adjust(ex_data)
#' ex_data = consensus(ex_data)
#' #write to a fasta file with the default output file
#' write_fastq(ex_data)
#' #specify the path to a custom output file, overwrite its current contents
#' write_fastq(ex_data, filename = 'example_out.fastq', append = FALSE)
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

  #TODO - need to make sure the phred scores are modified and carried through
  #then when the outseq is generated, also turn the phred numbers back into
  #the corresponding characters
  if(keep_phred == TRUE){
    if(is.null(x$phred)){
      stop("Cannot keep the phred scores for a DNAseq with no phred inputs")
    }
    outstring = paste(x$name, "\n",
                      x$outseq, "\n",
                      "+\n",
                      x$outphred, sep="")
  }else{
    outstring = paste(x$name, "\n",
                      x$outseq, "\n",
                      "+\n",
                      paste(rep(phred_placeholder, times = nchar(x$outseq)), collapse = ""), sep="")
  }
  
  write(outstring, file =  filename, append = append)
}

#' A wrapper function to deploy the fastq and fata output functions.
#' 
#' @param filename The name of the file to output the data to. Default is "denoised.fasta".
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

#' @rdname write_fastq
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
      write_fastq(x, ambig_char= ambig_char,
                  phred_placeholder = phred_placeholder, 
                  append = append, ...)
    }else{
      write_fastq(x, ambig_char = ambig_char,
                  filename = filename, 
                  phred_placeholder = phred_placeholder, 
                  append = append, ...)
    }
  }
  x
}