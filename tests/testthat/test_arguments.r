
test_that("Input arguments are handled correctly.", {
  
  #load in a set of example sequences  
  fastq_gz_real_file = system.file('extdata/coi_sequel_data_subset.fastq.gz', package = 'debar')
  data = read_fastq(fastq_gz_real_file)
  
  #
  ex1 = denoise(data$sequence[[1]], keep_phred = FALSE, to_file = FALSE)
  #x[['adjustment_count']] ==  1
  expect_equal(ex1[['adjustment_count']], 1)
  #is.null(ex1[['phred']]) == TRUE
  expect_equal(is.null(ex1[['phred']]), TRUE)
  
  super_short = 'GGTAGGAGCGTGTATACAGCGGCAGTCGAACATGTAGCTGACTCAGGTCACATTCAACAAATCATAAAGATAT'
  ex2 = denoise(super_short, keep_phred = FALSE, to_file = FALSE)
  
  #ex2$reject == TRUE
  expect_equal(ex2$reject, TRUE)
  
  ex3 = denoise(super_short, keep_phred = FALSE, terminate_rejects = FALSE, to_file = FALSE)
  #ex3$reject== FALSE
  expect_equal(ex3$reject, FALSE)
  #ex3$outseq == toupper(rev_comp(super_short))
  expect_equal(ex3$outseq, toupper(rev_comp(super_short)))
  
  #test that the warning is thrown if different lengths in the phred and DNA sequences.
  expect_error(denoise(data$sequence[[78]], phred = data$quality[[187]], to_file = FALSE), 
  "The length of the input DNA sequence and phred scores must match.")

})
