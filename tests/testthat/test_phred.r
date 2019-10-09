test_that("fastq phred scores can be processed and maintained correctly", {
  fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')
  
  test_data = read_fastq(fastq_gz_real_file, keep_quality = TRUE)

  test_data$phred = lapply(test_data$quality, seqinr::s2c)
  
  test_data$test_seq = lapply(test_data$sequence, seqinr::s2c)
  
  test_phred = test_data$phred[[1]]
  test_seq = test_data$test_seq[[1]]
  
  names(test_seq) = test_phred
  
  test_seq[1:5]
  
  #get the labelled bo
  test_seq[1]
  #get just the bp
  test_seq[[1]]
  #change the phred score for the bp
  names(test_seq)[1] = "*"
  #get the bp and the adjusted phred
  test_seq[1]
  
  test_building = test_seq[1:5]
  test_building = c(test_building, "N")
  names(test_building)[length(test_building)] = "*"
  names(test_building)=NULL
  
      
})