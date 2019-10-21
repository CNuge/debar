
test_that("Input arguments are handled correctly.", {
  
#load in a set of example sequences  
fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')
data = read_fastq(fastq_gz_real_file)

  
x = denoise(data$sequence[[71]], keep_phred = FALSE, to_file = FALSE)

#x[['adjustment_count']] ==  2
expect_equal(x[['adjustment_count']], 2)

#test that the warning is thrown if different lengths in the phred and DNA sequences.
expect_error(denoise(data$sequence[[78]], phred = data$quality[[187]], to_file = FALSE), 
"The length of the input DNA sequence and phred scores must match.")


})