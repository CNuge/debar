
#' Run the denoiser pipeline for a sequence read.
#'
#' 
#' @param x a DNA sequence string.
#' @param ... additional arguments to be passed between methods.
#' @param name an optional character string. Identifier for the sequence.
#' @param phred an optional character string. The phred score string corresponding to the nucleotide string.
#' If passed then the input phred scores will be modified along with the nucleotides and carried through
#' to the sequence output. Default = NULL.
#' @param phred_test Default = TRUE.
#' @param min_avg_qv
#' @param max_perc_low The maximum percentage of nucleotides in the string with QV values lower than 20. Default is 25%
#' @param max_perc_ultra_low The maximum percentage of nucleotides in the string with QV values lower than 10. Default is 5%
#' @param dir_check Should both the forward and reverse compliments be considered?
#' @param min_match The minimum number of sequential matches to the PHMM for a sequence to be denoised.
#' Otherwise flag the sequence as a reject.
#' @param censor_length the number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 5.
#' @param added_phred
#' @param added_phred The phred character to use for characters inserted into the original sequence.
#' @param adjust_limit the maximum number of corrections that can be applied to a sequence read. If this number is exceeded 
#' then the entire read is rejected. Default is 3.
#' @param keep_flanks Should the regions of the input sequence outside of the barcode region be readded to the denoised sequence
#' prior to outputting to the file. Default is TRUE. 
#' False will lead to only the denoised sequence for the 657bp barcode region being output to the file.
#' @param ambig_char The character to use for ambigious positions in the sequence that is output to the file. Default is N.
#' @param adjust_limit the maximum number of corrections that can be applied to a sequence read. If this number is exceeded 
#' then the entire read is masked with ambigious characters. Default is 5.
#' @param framed_output Boolean indicating if the output shouldhave leading dashes to establish a common
#' reading frame in the output sequences. Default is TRUE. Param only applied if keep_flanks = FALSE.
#' @param to_file Boolean indicating whether the sequence should be written to a file. Default is TRUE.
#' @param outformat The format of the output file. Options are fasta or fastq (default) format.
#' @param filename The name of the file to output the data to. Default filenames are respectively: denoised.fasta or denoised.fastq.
#' @param append Should the denoised sequence be appended to the output file?(TRUE) 
#' Or should the sequence overwrite the output file?(FALSE) Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'. 
#' Used with write_fastq and keep_phred == FALSE only.
#' @param terminate_rejects Boolean indicating if analysis of sequences that fail to meet phred quality score or path 
#' match thresholds should be terminated early (prior to sequence adjustment and writing to file). Default it true.
#' @param aa_check Boolean indicating whether the amino acid sequence should be generated and assessed for stop codons. Default = FALSE.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids. Default is 0, meaning
#' that censored translation is performed (amigious codons ignored). Used only when aa_check = TRUE.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would offset the sequence by one and therefore make the second character in the
#' framed sequence the the first nucelotide of the first codon. Used only when aa_check = TRUE.
#' @return a class object of code{"DNAseq"} 
#' @examples
#' # Denoise example sequence with default paramaters.
#' ex_data = denoise(example_nt_string, name = 'example_sequence_1')
#' @export
#' @name denoise
denoise = function(x, ...){
  UseMethod("denoise")
} 

#' @rdname denoise
#' @export
denoise.default = function(x, ...,
                             name = character(),
                             phred = NULL, 
                             phred_check = TRUE,
                             min_avg_qv = 20,
                             max_perc_low = 0.25,
                             max_perc_ultra_low = 0.05,
                             dir_check = TRUE, 
                             min_match = 100,
                             censor_length = 5,
                             added_phred = "*",
                             adjust_limit = 5,
                             ambig_char = "N",
                             to_file = TRUE,
                             keep_flanks = TRUE,
                             keep_phred = TRUE,
                             outformat = "fastq",
                             terminate_rejects = TRUE,
                             filename = NULL, 
                             phred_placeholder = "#",
                             write_masked = TRUE,
                             aa_check = FALSE, 
                             trans_table = 0,
                             frame_offset = 0,
                             append = TRUE
                             ){

  dat = DNAseq(x, name = name , phred = phred)

  if(phred_test==TRUE){
    dat = phred_check(dat, min_avg_qv = min_avg_qv, 
                           max_perc_low = max_perc_low, 
                           max_perc_ultra_low=max_perc_ultra_low, 
                           ...)  
  }else{
    dat$reject = FALSE
  }
  
  if(dat$reject == TRUE && terminate_rejects == TRUE){
    return(dat)
  }

  dat = frame(dat, dir_check = dir_check, 
                   min_match = min_match,
                   ...)
  
  if(dat$reject == TRUE && terminate_rejects == TRUE){
    return(dat)
  }
  
  dat = adjust(dat, censor_length= censor_length,
                    added_phred = added_phred,
                    ...)
  if(aa_check == TRUE){
    dat = aa_check(dat, trans_table = trans_table, 
                        frame_offset = frame_offset,
                        ...)
  }
  
  dat = outseq(dat, keep_flanks = keep_flanks, 
                    ambig_char = ambig_char, 
                    adjust_limit = ambig_char,
                    ...)
  
  if(dat$reject == TRUE && terminate_rejects == TRUE){
    return(dat)
  }
  
  if(to_file == TRUE){
    dat = write_wrapper(dat, filename = filename, 
                             outformat = outformat,
                             append=append,
                             keep_phred = keep_phred,
                             phred_placeholder = phred_placeholder,
                             ...)
  } 
  dat
}


#' Run the log file and reject processing code after running the denoise pipeline.
#'
#' @keywords internal
meta_check = function(temp, log_data = list(), keep_rejects = FALSE, log_file = FALSE, ...){
  #write read to the reject file if keep_rejects option enabled
  if(keep_rejects == TRUE && temp$reject == TRUE){
    temp$outseq = temp$raw
    temp$outphred = temp$phred
    write_wrapper(temp, ...)
  }
  #add to the summary stats if the log file option is enabled
  if(log_file == TRUE){
    log_data[['total_reads']] = log_data[['total_reads']] + 1
    if(x$reject == TRUE){
      log_data[['reject_count']] = log_data[['reject_count']] + 1
    }else{
      log_data[['good_count']] = log_data[['good_count']] + 1
    }
    if(x$adjustment_count>0){
      log_data[['good_denoised']] = log_data[['good_denoised']] + 1
    }else{
      log_data[['good_unaltered']] = log_data[['good_unaltered']] + 1 
    }
  }
  log_data
}


#'Denoise the sequence data from a given file.
#'
#'@param x The name of the file to denoise sequences from.
#'@param filename The name of the file to output sequences to.
#'@param file_type The format of the file to be denoised. Options are fastq or fasta. Default is fastq.
#'@param log_file Boolean indicating if a log file should be produced. Default is FALSE.
#'@param append_log, Boolean indicating if a log data should be appended to an existing log file. Default is TRUE.
#'@param keep_rejects Boolean indicating if the bad reads should be written to a separate file (with the name
#'"rejects_" + filename). Defaut is false
#'@param multicore An integer specifying the number of cores over which to multithread the denoising process. 
#'Default is FALSE, meaning the process is not multithreaded.
#'@param ... additional arguments to be passed to the \link{denoise} function.
#'
#'@seealso \code{\link{denoise}}
#'@examples
#' fastq_dat_file = system.file('extdata/ccs_subset.fastq', package = 'coiDenoiser')
#' denoise_file(fastq_dat_file, to_file = FALSE)
#'
#' @export
#' @name denoise_file  
denoise_file = function(x, ...){
    UseMethod("denoise_file")
  } 

#' @rdname denoise_file
#' @export
denoise_file.default = function(x, ..., filename = 'output.fastq',  file_type = "fastq", 
                                  log_file = FALSE, keep_rejects = FALSE, multicore = FALSE){
  #set up additional output paramaters if needed
  if(keep_rejects == TRUE){
    reject_filename = paste0("rejects_",filename)
  }
  if(log_file == TRUE){
    log_data = list(input_file = x, total_reads = 0, good_count = 0, reject_count = 0, good_denoised = 0, good_unaltered = 0,  
                    start_time = Sys.time(), end_time = NA, time_elapsed = NA)
    
    log_filename = paste0("log_",filename)  
  }
  
  #read in the data
  if(file_type == "fastq"){
    data = read_fastq(x)
  }else if(file_type == "fasta"){
    data = read_fasta(x)
  }else{
    stop("file_type must be either fasta or fastq")
  }

  if(multicore == FALSE){
    for(i in 1:length(data$sequence)){
      temp = denoise(data$sequence[[i]], filename = filename, name = data$header_data[[i]], phred = data$quality[[i]], ...)
      log_data = meta_check(temp, ...)
    }
  }else{
    parallel::mclapply(1:length(data$sequence), function(i){
      temp = denoise(data$sequence[[i]], filename = filename, name = data$header_data[[i]], phred = data$quality[[i]], ...)
      #TODO - if the reassignment from within the mclapply leads to errors
      #then remove the log_data assignment here. still need the bad reads saved
      log_data = meta_check(temp, log_data, ...)}, mc.cores = multicore)
  }
  
  if(log_file == TRUE){
    log_data[['end_time']] = Sys.time()
    log_data[['time_elapsed']] = difftime(log_data[['end_time']] , log_data[['start_time']] , units = "mins")
    log_data= data.frame(log_data)
    write.csv(log_data, log_filename, row.names = FALSE, row.names=FALSE, append = append_log)
  }
}


