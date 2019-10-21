test_that("Sequence files are read in properly.", {
  
  fastq_dat_file = system.file('extdata/ccs_subset.fastq', package = 'seqdenoise')
  data1 = read_fastq(fastq_dat_file)
  
  #all.equal(dim(data1), c(50,3))
  expect_equal(all.equal(dim(data1), c(50,3)), TRUE)
  
  fastq_real_file = system.file('extdata/sequel_smrt_subset.fastq', package = 'seqdenoise')
  data2 = read_fastq(fastq_real_file)
  
  #all.equal(dim(data2), c(250,3))
  expect_equal(all.equal(dim(data2), c(250,3)), TRUE)
  
  fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')
  data3 = read_fastq(fastq_gz_real_file, keep_quality = FALSE)
  
  #all.equal(dim(data3), c(250,2))
  expect_equal(all.equal(dim(data3), c(250,2)), TRUE)
  
})