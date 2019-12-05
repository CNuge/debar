test_that("Amino acid sequences are checked properly.", {

  ex_data = DNAseq(example_nt_string, name = 'ex1')
  ex_data =  frame(ex_data)
  ex_data = adjust(ex_data)
  
  #test with no translation table
  ex_data1 = aa_check(ex_data)
  expect_equal(ex_data1$reject, FALSE)
  
  #test with  translation table
  ex_data2 = aa_check(ex_data, trans_table = 2)
  expect_equal(ex_data2$reject, FALSE)
  
  #use the wrong translation table to force an aa_check flag
  ex_data3 = aa_check(ex_data, trans_table = 11)
  expect_equal(ex_data3$reject, TRUE)
  
  #cleanup tests for the hard to hit lines
  expect_equal(translate_codon("NNA"), "-")
  expect_equal(translate_codon("Ann"), "-")
  
  })