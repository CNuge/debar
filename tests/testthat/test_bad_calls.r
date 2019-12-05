test_that("The correct errors are thrown for bad denoise calls.", {
  
  #non DNA characters passed to the DNAseq object
  err_call = "Unallowed character in DNA string: d 
Valid characters are: a t g c - n"
  expect_error(DNAseq("IamnotDNA"), err_call)

  expect_error(DNAseq("ATGC", name = "mismatch_size", phred = "~~~~~~~~~~"), 
               "The length of the input DNA sequence and phred scores must match.")

  expect_error(DNAseq(),  "Must pass a DNA sequence.")    
})