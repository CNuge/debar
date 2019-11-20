test_that("Sequence files are read in properly.", {
  
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
  
})