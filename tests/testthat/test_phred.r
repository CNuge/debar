test_that("fastq phred scores can be processed and maintained correctly", {
  fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')
  
  test_data = read_fastq(fastq_gz_real_file, keep_quality = TRUE)

  test_data$num_phred = lapply(test_data$quality, phred2numeric)
  
  test_phred = test_data$num_phred[[1]]
  test_phred = unname(test_phred)
  
  
  test_seq = seqinr::s2c(test_data$sequence[[1]])
  
  names(test_seq) = test_phred
  
  test_seq[1:5]
  
  test_seq[[1]]
  
})