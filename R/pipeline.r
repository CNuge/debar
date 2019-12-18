
#' Run the denoiser pipeline for a sequence read.
#'
#' This function runs the complete denoising pipeline for a given input sequence and its corresponding
#' name and phred scores. The default behaviour is set to interface with fastq files (standard output for
#' most sequencers). 
#' 
#' Since the pipeline is designed for recieving or outputting either fasta or fastq data, this function is 
#' hevaily paramaterized. Note that not all paramaters will affect all use cases (i.e. if your outformat is to 
#' a fasta file, then the phred_placeholder paramater is ignored
#' as this option only pertains to fastq outputs). The user is encouraged to read the vignette for a detailed 
#' walkthrough of the denoiser pipeline that will help identify the paramaters that relate to their given needs.
#' 
#' @param x a DNA sequence string.
#' @param ... additional arguments to be passed between methods.
#' @param name an optional character string. Identifier for the sequence.
#' @param phred an optional character string. The phred score string corresponding to the nucleotide string.
#' If passed then the input phred scores will be modified along with the nucleotides and carried through
#' to the sequence output. Default = NULL.
#' @param keep_phred Should the original PHRED scores be kept in the output? Default is TRUE.
#' @param dir_check A boolean indicating if both the forward and reverse compliments of a sequence should 
#' be checked against the PHMM. Default is TRUE.
#' @param double_pass A boolean indicating if a second pass through the Viterbi algorithm should be conducted for sequences
#' that had leading nucleotides not matching the PHMM. This improves the accurate establishment of reading frame and will
#' reduce false rejections by the amino acid check, but this comes at a cost of additional processing time. Default is TRUE.
#' @param min_match The minimum number of sequential matches to the PHMM for a sequence to be denoised.
#' Otherwise flag the sequence as a reject.
#' @param max_inserts The maximum number of sequention insert states occuring in a sequence 
#' (including the flanking regions). If this number is
#' exceeded than the entire read will be discarded if terminate_rejects = TRUE. Default is 400.
#' @param censor_length the number of base pairs in either direction of a PHMM correction
#' to convert to placeholder characters. Default is 7.
#' @param added_phred The phred character to use for characters inserted into the original sequence.
#' @param adjust_limit the maximum number of corrections that can be applied to a sequence read. If this number is exceeded 
#' then the entire read is rejected. Default is 3.
#' @param keep_flanks Should the regions of the input sequence outside of the barcode region be readded to the denoised sequence
#' prior to outputting to the file. Options are TRUE, FALSE and 'right'. The 'right' option will keep the trailing flank
#' but remove the leading flank. Default is TRUE. 
#' False will lead to only the denoised sequence for the 657bp barcode region being output to the file.
#' @param ambig_char The character to use for ambigious positions in the sequence that is output to the file. Default is N.
#' @param to_file Boolean indicating whether the sequence should be written to a file. Default is TRUE.
#' @param outformat The format of the output file. Options are fasta or fastq (default) format.
#' @param outfile The name of the file to output the data to. Default filenames are respectively: denoised.fasta or denoised.fastq.
#' @param append Should the denoised sequence be appended to the output file?(TRUE) 
#' Or should the sequence overwrite the output file?(FALSE) Default is TRUE.
#' @param phred_placeholder The character to input for the phred score line. Default is '#'. 
#' Used with write_fastq and keep_phred == FALSE only.
#' @param terminate_rejects Boolean indicating if analysis of sequences that fail to meet phred quality score or path 
#' match thresholds should be terminated early (prior to sequence adjustment and writing to file). Default it true.
#' @param aa_check Boolean indicating whether the amino acid sequence should be generated and assessed for stop codons. Default = TRUE.
#' @param trans_table The translation table to use for translating from nucleotides to amino acids. Default is 0, meaning
#' that censored translation is performed (amigious codons ignored). Used only when aa_check = TRUE.
#' @param frame_offset The offset to the reading frame to be applied for translation. By default the offset
#' is zero, so the first character in the framed sequence is considered the first nucelotide of the first codon.
#' Passing frame_offset = 1 would offset the sequence by one and therefore make the second character in the
#' framed sequence the the first nucelotide of the first codon. Used only when aa_check = TRUE.
#' @return a class object of code{"DNAseq"} 
#' @examples
#' # Denoise example sequence with default paramaters.
#' ex_data = denoise(example_nt_string_errors, 
#'                   name = 'example_sequence_1', 
#'                   keep_phred = FALSE, 
#'                   to_file = FALSE)
#' 
#' #fastq data from a file
#' #previously run
#' fastq_example_file = system.file('extdata/coi_sequel_data_subset.fastq', 
#'                                  package = 'debar')
#' data = read_fastq(fastq_example_file)
#' #denoise the first sequence in the file
#' #use a custom censor length and no amino acid check
#' dn_dat_1 = denoise(x = data$sequence[[1]], 
#'                     name = data$header[[1]], 
#'                     phred = data$quality[[1]], 
#'                     censor_length = 11, 
#'                     aa_check = FALSE, 
#'                     to_file = FALSE)
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
                             dir_check = TRUE,
                             double_pass = TRUE,
                             min_match = 100,
                             max_inserts = 400,
                             censor_length = 7,
                             added_phred = "*",
                             adjust_limit = 5,
                             ambig_char = "N",
                             to_file = FALSE,
                             keep_flanks = TRUE,
                             keep_phred = TRUE,
                             outformat = "fastq",
                             terminate_rejects = TRUE,
                             outfile = NULL, 
                             phred_placeholder = "#",
                             aa_check = TRUE, 
                             trans_table = 0,
                             frame_offset = 0,
                             append = TRUE
                             ){

  dat = DNAseq(x, name = name , phred = phred)
  if(is.null(dat$phred)){
    keep_phred = FALSE
  } 
  
  dat = frame(dat, dir_check = dir_check,
                   double_pass = double_pass,
                   min_match = min_match,
                   max_inserts = max_inserts,
                   terminate_rejects= terminate_rejects,
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

  if(dat$reject == TRUE && terminate_rejects == TRUE){
    return(dat)
  }

  dat = outseq(dat, keep_flanks = keep_flanks, 
                    ambig_char = ambig_char, 
                    adjust_limit = ambig_char,
                    ...)
  
  if(dat$reject == TRUE && terminate_rejects == TRUE){
    return(dat)
  }
  
  if(to_file == TRUE){
    dat = write_wrapper(dat, outfile = outfile, 
                             outformat = outformat,
                             append = append,
                             keep_phred = keep_phred,
                             phred_placeholder = phred_placeholder,
                             ...)
  } 
  dat
}


#' Run the log file and reject processing code after running the denoise pipeline.
#'
#' @keywords internal
meta_check = function(x, log_data = list(), log_file = FALSE, keep_rejects = FALSE, reject_filename = "rejects.fastq", ...){
  #write read to the reject file if keep_rejects option enabled
  if(keep_rejects == TRUE && x$reject == TRUE){
    x$outseq = x$raw
    x$outphred = x$phred
    write_wrapper(x, outfile = reject_filename, ...)
  }
  #add to the summary stats if the log file option is enabled
  if(log_file == TRUE){
    log_data[['total_reads']] = log_data[['total_reads']] + 1
    if(x$reject == TRUE){
      log_data[['reject_count']] = log_data[['reject_count']] + 1
    }else{
      log_data[['good_count']] = log_data[['good_count']] + 1
      
      if(x$adjustment_count>0){
        log_data[['good_denoised']] = log_data[['good_denoised']] + 1
      }else{
        log_data[['good_unaltered']] = log_data[['good_unaltered']] + 1 
      }
    }
  }
  log_data
}


#' Denoise sequence data from a given file.
#'
#' This function allows for direct input to output exectution of the denoising pipeline. All paramaters for the 
#' fasta/fastq input and output functions as well as the denoise pipeline can be passed to this function. Please consult
#' the documentation for those functions for a list of available paramaters. The function will
#' run the denoise pipeline with the specified paramaters for all sequences in the input file, and write the denoised sequences
#' and corresponding header/quality information to the output file. Additionally the function allows for rejected reads to
#' be kept and sequestered to an additional output file (as opposed to being discarded) and also allows for a log file to
#' be produced that tracks several statistics including the execition time, number of denoised reads and number of
#' rejected reads.
#' 
#' Using this function is optimized by the appropriation of the multicore option, which allows a user to specify a number of
#' cores that the denoising process should be multithreaded across. The more cores available, the faster the denoising of the 
#' input data. It should be noted that the multithreading relies on the entire fastq file being read into memory, because of
#' this your machine's available ram will need to exceed the size of the unzipped fastq file being denoised. If your file
#' size exceeds the available memory you may want to consider spliting the input into several smaller files and denoising them
#' each with this function (this is a fast solution as the multicore option can be used to speed up denoising). Alternatively,
#' you can depoly the `denoise` function in an iterative fashion, reading in/denoising and writing only a single fq entry at
#' a time. This would require a much smaller memory footprint, but would be much slower due to the lack of multithreading.
#' 
#' @param x The name of the file to denoise sequences from.
#' @param outfile The name of the file to output sequences to.
#' @param informat The format of the file to be denoised. Options are fastq or fasta. Default is fastq.
#' @param outformat The format of the output file. Options are fasta or fastq (default) format.
#' @param to_file Boolean indicating whether the sequence should be written to a file. Default is TRUE.
#' @param log_file Boolean indicating if a log file should be produced. Default is FALSE.
#' @param keep_rejects Boolean indicating if the bad reads should be written to a separate file (with the name
#' "rejects_" + outfile). Defaut is FALSE.
#' @param multicore An integer specifying the number of cores over which to multithread the denoising process. 
#' Default is FALSE, meaning the process is not multithreaded.
#' @param ... additional arguments to be passed to the \link{denoise} and input/output functions.
#'
#' @seealso \code{\link{denoise}}
#'
#' @export
#' @name denoise_file  
denoise_file = function(x, ...){
    UseMethod("denoise_file")
  } 

#' @rdname denoise_file
#' @export
denoise_file.default = function(x, ..., outfile = 'output.fastq',  informat = "fastq", outformat = 'fastq',
                                  to_file = TRUE, log_file = FALSE, keep_rejects = FALSE, multicore = FALSE){
  #set up additional output paramaters if needed
  if(keep_rejects == TRUE){
    reject_filename = paste(c(strsplit(outfile,"\\.")[[1]][[1]], "_rejects.", strsplit(outfile,"\\.")[[1]][[2]]), collapse = "")  
  }
  if(log_file == TRUE){
    log_data = list(input_file = x, total_reads = 0, good_count = 0, reject_count = 0, good_denoised = 0, good_unaltered = 0,  
                    start_time = Sys.time(), end_time = NA, time_elapsed = NA)
    
    log_filename = paste(c(strsplit(outfile,"\\.")[[1]], "_log.out"), collapse = "")  
  }else{
    log_data = list()
  }
  
  #read in the data
  if(informat == "fastq"){
    print(paste0("Reading fastq file:", x))
    data = read_fastq(x)
  }else if(informat == "fasta"){
    print(paste0("Reading fasta file:", x))
    data = read_fasta(x)
  }else{
    stop("informat must be either fasta or fastq")
  }
  
  if(substr(outfile, nchar(outfile)-2, nchar(outfile)) == ".gz"){
    outfile=substr(outfile, 1, nchar(outfile)-3)
  }

  print(paste0("Denoising data from file"))
  if(multicore == FALSE){
    for(i in 1:length(data$sequence)){
      temp = denoise(data$sequence[[i]], to_file = to_file, 
                                          outfile = outfile,
                                          outformat = outformat,
                                          name = data$header_data[[i]], 
                                          phred = data$quality[[i]], ...)
      log_data = meta_check(x = temp, log_data = log_data, 
                                  keep_rejects = keep_rejects, 
                                  log_file = log_file,
                                  reject_filename = reject_filename, ...)
    }
  }else{
    print(paste0("multithreading across ", multicore, " cores."))
    if(log_file == FALSE){
      parallel::mclapply(1:length(data$sequence), function(i, ...){
        temp = denoise(data$sequence[[i]], to_file = to_file, 
                                            outfile = outfile,
                                            outformat = outformat,
                                            name = data$header_data[[i]], 
                                            phred = data$quality[[i]], ...)
        meta_check(x = temp, log_data = FALSE, 
                             keep_rejects = keep_rejects, 
                             log_file = log_file,
                             reject_filename = reject_filename, ...)
        NULL
      }, mc.cores = multicore)
    }else{
      log_row = list(total_reads = 0, good_count = 0, reject_count = 0, good_denoised = 0, good_unaltered = 0)
      log_rows = parallel::mclapply(1:length(data$sequence), function(i, ...){
        temp = denoise(data$sequence[[i]], outfile = outfile, outformat = outformat, name = data$header_data[[i]], phred = data$quality[[i]], ...)
        outrow = meta_check(x = temp, log_data = log_row, 
                   keep_rejects = keep_rejects, 
                   log_file = log_file,
                   reject_filename = reject_filename,  ...)
        data.frame(outrow)
      }, mc.cores = multicore)
      
      mc_log = do.call(rbind,log_rows)
      
      log_data[['total_reads']] = sum(mc_log$total_reads)
      log_data[['good_count']] = sum(mc_log$good_count)
      log_data[['reject_count']] = sum(mc_log$reject_count)
      log_data[['good_denoised']] = sum(mc_log$good_denoised)
      log_data[['good_unaltered']] = sum(mc_log$good_unaltered)

    }
  }
  
  if(log_file == TRUE){
    log_data[['end_time']] = Sys.time()
    log_data[['time_elapsed']] = difftime(log_data[['end_time']] , log_data[['start_time']] , units = "mins")
    log_data = data.frame(log_data)
    utils::write.csv(log_data, log_filename, row.names = FALSE)
  }
  
  print("Done.")
}



#' List-to-list denoising of COI barcode sequences.
#' 
#' This function provides a shortcut for running the denoise function
#' on a list of sequences. The to_return option can be used to control
#' whether this function returns a list of sequence strings (default), 
#' or a list of DNA seq objects.
#'
#'
#' @param x A list like object of barcode sequences.
#' @param to_return Indicate whether a the function should return a list of 
#' sequence ('seq') or the full DNAseq object ('DNAseq). Default is ('seq')
#' @param cores The number of cores across which to thread the denosiing. Default is 1.
#' @param ... additional arguments to pass to the denoise algorithm.
#' @seealso \code{\link{denoise}}
#' @examples
#' #denoise a list of sequences
#' out = denoise_list(ex_nt_list, dir_check = FALSE, double_pass = FALSE)
#' #denoise and add placehers to outputs 
#' 
#' #return a list of DNAseq objects 
#' ex_DNAseq_out = denoise_list(ex_nt_list, to_return = 'DNAseq',
#'  dir_check = FALSE, double_pass = FALSE)
#'
#' @export
#' @name denoise_list  
denoise_list = function(x, to_return = 'seq', cores = 1, ...){
  
  out = parallel::mclapply(x, function(y){
    denoise(y, to_file = FALSE, keep_phred = FALSE, ...)
  }, mc.cores = cores)
  
  if(to_return == 'seq'){
    seqs = lapply(out, function(y){
      y[['outseq']]
    })
    
    return(seqs)
  }
  
  return(out)
  
}
