test_that("The file to file denoising pipeline performs as expected.", {
  #test the full pipeline on a small sequence datasets
  # have the 
  fastq_test_file = system.file('extdata/small_unittest.fastq', package = 'debar')
  fasta_test_file = system.file('extdata/small_unittest.fasta', package = 'debar')
  
  #
  # Write to a fastq file
  #
  temp3 = file.path(tempdir(), "temp_out.fastq")
  on.exit(unlink(temp3))
  
  
  denoise_file(fastq_test_file, outfile = temp3,  
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)

  #test multicore
  denoise_file(fastq_test_file, outfile = temp3,  multicore = 2,
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  denoise_file(fasta_test_file, outfile = temp3,  multicore = 2,
               informat = 'fasta', outformat="fastq",
               keep_phred = FALSE, phred_placeholder = "~",
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  denoise_file(fasta_test_file, outfile = temp3, informat = 'fasta', outformat="fastq",
               keep_phred = FALSE, phred_placeholder = "~", 
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  #
  # Write to a fasta file
  #
  temp4 = file.path(tempdir(), "temp_out.fasta")  
  temp_rej = file.path(tempdir(), "temp_out_rejects.fasta")
  temp_log = file.path(tempdir(), "temp_out_log.out")
  on.exit(unlink(temp4))
  on.exit(unlink(temp_rej))
  on.exit(unlink(temp_log))
  
  denoise_file(fasta_test_file, outfile = temp4, informat = 'fasta', outformat="fasta",
               keep_phred = FALSE, phred_placeholder = "~", 
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  denoise_file(fastq_test_file, outfile = temp4, informat = "fastq",
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  temp4 = 'out.fq'
  
  #other param combos
  denoise_file(fastq_test_file, outfile = temp4, informat = "fastq",
               log_file = FALSE, keep_rejects = TRUE, double_pass = FALSE, dir_check = FALSE)

  denoise_file(fastq_test_file, outfile = temp4, informat = "fastq",
               log_file = TRUE, keep_rejects = TRUE, double_pass = FALSE, dir_check = FALSE)
  
  denoise_file(fastq_test_file, outfile = temp4, informat = "fastq",  multicore = 2,
               log_file = TRUE, keep_rejects = TRUE, double_pass = FALSE, dir_check = FALSE)
  
  # Test bad calls of file types
  expect_error(denoise_file(fastq_test_file, outfile = temp4, informat="text",
                            log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE),
               "informat must be either fasta or fastq")
  
  expect_error(denoise_file(fasta_test_file, outfile = temp4, informat = 'fasta', outformat="text",
                            keep_phred = FALSE, phred_placeholder = "~", 
                            log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE),
               "Invalid output format! Must be one of: 'fasta', 'fastq' or 'none'")
  
  
})