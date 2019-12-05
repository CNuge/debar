test_that("Sequence files are read, denoised and written properly.", {
  
  fasta_real_file = system.file('extdata/coi_sequel_data_subset.fasta', package = 'debar')
  data1 = read_fasta(fasta_real_file)
  
  #all.equal(dim(data1), c(50,2))
  expect_equal(all.equal(dim(data1), c(250,2)), TRUE)
  
  fastq_real_file = system.file('extdata/coi_sequel_data_subset.fastq', package = 'debar')
  data2 = read_fastq(fastq_real_file)
  
  #all.equal(dim(data2), c(250,3))
  expect_equal(all.equal(dim(data2), c(250,3)), TRUE)
  
  fastq_gz_real_file = system.file('extdata/coi_sequel_data_subset.fastq.gz', package = 'debar')
  data3 = read_fastq(fastq_gz_real_file, keep_quality = FALSE)
  
  #all.equal(dim(data3), c(250,2))
  expect_equal(all.equal(dim(data3), c(250,2)), TRUE)

  #test writing to fastq file using a temp
  temp = tempfile()
  on.exit(unlink(temp))
  
  #write to a fastq when no quality data
  x1 = denoise(data1$sequence[[1]], data1$header[[1]], double_pass = FALSE, dir_check = FALSE)
  write_fastq(x1 , outfile = temp, keep_phred = FALSE, phred_placeholder = "~")
  
  #write to fastq with quality data
  x2 = denoise(data2$sequence[[1]], name = data2$header[[1]], phred = data2$quality[[1]], double_pass = FALSE, dir_check = FALSE)
  write_fastq(x2 , outfile = temp)
  
  #test writing to fasta using a temp
  temp2 = tempfile()
  on.exit(unlink(temp2))
  
  #write to a fastq when no quality data
  x3 = denoise(data1$sequence[[1]], data1$header[[1]], double_pass = FALSE, dir_check = FALSE)
  write_fastq(x3 , outfile = temp2, keep_phred = FALSE, phred_placeholder = "~")
  
  #write to fastq with quality data
  x4 = denoise(data2$sequence[[1]], name = data2$header[[1]], phred = data2$quality[[1]], double_pass = FALSE, dir_check = FALSE)
  write_fastq(x4 , outfile = temp2)
  
  #test the full pipeline on a small sequence datasets
  # have the 
  fastq_test_file = system.file('extdata/small_unittest.fastq', package = 'debar')
  fasta_test_file = system.file('extdata/small_unittest.fasta', package = 'debar')
  
  #
  # Write to a fasta file
  #
  temp3 = file.path(tempdir(), "temp_out.fastq")
  on.exit(unlink(temp3))
  
  denoise_file(fastq_test_file, outfile = temp3,  
                log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  denoise_file(fasta_test_file, outfile = temp3, file_type = 'fasta', keep_phred = FALSE, phred_placeholder = "~", 
              log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
  
  #
  # Write to a fasta file
  #
  temp4 = file.path(tempdir(), "temp_out.fasta")  
  on.exit(unlink(temp4))
  
  denoise_file(fasta_test_file, outfile = temp4, file_type = 'fasta', keep_phred = FALSE, phred_placeholder = "~", 
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)

  denoise_file(fastq_test_file, outfile = temp4, 
               log_file = FALSE, keep_rejects = FALSE, double_pass = FALSE, dir_check = FALSE)
  
})
