
test_that("Input arguments are handled correctly.", {
  
#load in a set of example sequences  
fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')
data = read_fastq(fastq_gz_real_file)

  
x = denoise(data$sequence[[71]], keep_phred = FALSE, to_file = FALSE)
x



x = denoise(data$sequence[[78]], phred = data$quality[[187]])
x


})