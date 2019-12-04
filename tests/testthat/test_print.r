test_that("DNAseq objects are printed correctly.", {
  
  toy_example1 = 'ACAATGTATTTCATCCTAGGAATATGATCAGGAATATTAGGAATAATATTAAGAATATTATTCGAATTGAACTAGCACAACCAGGACCATTAATAAGAAATGACCAAATTTATAATGTAATTGTTACATCTCATGCATTCATTATAATTTTTATGGTAATACCAATTATAATCGGAGGATTCGGAAATTGATTAGTACCATTAATAATTGGAGCACCAGATATAGCATTCCCACGAATAAATAATATAAGATTTTGACTTCTTCCACCCTCACTTACATTGTTAATTTCAAGATCAATAGTAGAAATAGGACCAGGAACAGGATGAACACTATATCCACCGCTATCGTCAAATATTGCACATTCGGGAGGAAGTGTAGACTTAAGAATTTTTTCATTACACTTAGCCGGTGTATCATCAATCCTAGGAGCAATTAACTTTATTACAACAATTCTAAATATACGAACACTAGGAATATCATTAGACCGAACCCCCTTATTCGTATGATCAGTAATAATTACTGCAATTTTACTACTTCTATCTCTACCAGTACTAGCTGGAGCAATCACGATA'
  test_toy1_out = denoise(toy_example1, 
                        name = "example1",
                        ambig_char = "N",
                        censor_length = 3,
                        dir_check = FALSE,
                        aa_check = TRUE, 
                        to_file = FALSE)


  expected = "A DNAseq object.
Raw Sequence:
acaatgtatttcatcctaggaatat...agtactagctggagcaatcacgata
Corrections applied: 1
Output Sequence:
ACAATGTATTTCATCCTAGGAATAT...AGTACTAGCTGGAGCAATCACGATA"

  outstring = capture_output(test_toy1_out, print=TRUE)

  expect_equal(outstring, expected)
})