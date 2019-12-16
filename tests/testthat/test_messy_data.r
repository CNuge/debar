test_that("Messy SEQUEL sequences are processed as anticipated", {
  
  ############################################################################
  # Real read test 1
  ############################################################################
  #this has a strong run of barcode data, but appears to be chimeric/ a half fwd half back
  #this one should be rejected
  
  t1 = "gctcctaaaattgaagaaactcctgctaaatgtaaagaaaaaattgctaaatctactgaagcaccaccatgagcaataatagaagataaaggtggataaacagttcatccagtaccagccccattttctactatacttctagccattaatagagtcaatgaaggcggtaataatcagaaactcatattatttattcgaggaaatgctatatccggagctcctaatattaatggaactaatcaatttccaaatcctccaattataattggtattactataaaaaaaattattacgaaagcatgggcagtgacaattacattataaatttgatcatcaccaattaatgctcctggatgtcctaattcagctcgaattaagattcttagagaagttcctactattccagctcatgctccaaatataaagtacaaagttccaatatctttatgattggttgactgtgacctgagtcagctacatgttcgactgcatagacagtgcgctgtctatgcagtcgaacatgtagctgactcaggtcacggtcaacaaatcataaagatattggaactttgtactttatatttggagcatgagctggaatagtaggaacttctctaagaatcttaattcgagctgaattaggacatccaggagcattaattggtgatgatcaaatttataatgtaattgtcactgcccatgctttcgtaataatttttttatagtaataccaattataattggaggatttggaaattgattagttccattaatattaggagctccagatatagcatttcctcgaataaataatatgagtctctgattattaccgccttcattgactctattaatggctagaagtatagtagaaaatggggctggtactggatgaactgtttatccacctttatcttctattattgctcatggtggtgcttcagtagatttagcaattttttctttacatttagcaggagtttcttcaattttaggagc"
  
  test_seqdenoise = DNAseq(t1, name = "test1")
  test_seqdenoise = frame(test_seqdenoise)
  #test_seqdenoise$data$ntPHMMout[['path']]
  #test_seqdenoise$reject == TRUE
  expect_equal(test_seqdenoise$reject, TRUE)
  
  
  ############################################################################
  # Real read test 2
  ############################################################################
  
  #issue here, the 296 characters outside of the denoised area are dropped
  #i.e. after the run of 2s happens, all of the following data is omitted in the subsequent output.
  #need adjust seq to keep these parts on the edges, find the logic that ceases the iteration if 
  # a large string of 1s or twos is encountered... consider how to proceed given that the later
  #example, t3 is benefitting from the trunction at this early point.
  t2 = "CGTATGTATGTCGCGCAGTCGAACATGTAGCTGACTCAGGTCACGGTCAACAAATCATAAAGATATTGGAACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCTTTAAGTTTATTAATTCGACTTGAATTAAGAACAATAAAAAATTTAATTGGTAATGATCAAATTTATAATGTAATTGTTACAGCCCACGCATTTATTATAATTTTTTTTATAGTTATACCAATTTTAATTGGAGGATTTGGAAATTGATTAGTCCCTATTATATTAGGAGCCCCAGATATAGCTTTTCCACGAATAAATAATATAAGATTCTGAATATTACCCCCTTCTTTATCTTTACTATCAATTAGAAGAATAGTAGAAACTGGGACCGGAACAGGATGAACTGTTTACCCACCTCTTTCTTCTGTAATTGCACATACTGGATCTTCAGTAGATTTTTCAATTTTTTCTCTACATATTGCAGGAATTTCATCTATTTTAGGAGCAATTAATTTTATTTCTACAATAATAAATATAAAAATTGAAAAATCTACTGAAGATCCAGTATGTGCAATTACAGAAGAAAGAGGTGGGTAAACAGTTCATCCTGTTCCGGTCCCAGTTTCTACTATTCTTCTAATTAATAGTAAAGATAAAGAAGGGGGTAATATTCAGAATCTTATATTATTTATTCGTGGAAAAGCTATATCTGGGGCTCCTAATATAATAGGGACTAATCAATTTCCAAATCCTCCAATTAAAATTGGTATAACTATAAAAAAAAATTATAATAAATGCGTGGGCTGTAACAATTACATTATAAATTTGATCATTACCAATTAAATTTTTTATTGTTCTTAATTCAAGTCGAATTAATAAACTTAAAGAAGTTCCTACTATTCCTGCTCAAATTCCAAAAATAAAATATAAAGTTCCAATATCTTTATGATTGGTTGAATGTGACCTGAGTCAGCTACATGTTCGACTGCGCGACATACATACG"
  
  test_seqdenoise2 = DNAseq(t2, name = "test1")
  test_seqdenoise2 = frame(test_seqdenoise2)
  test_seqdenoise2 = adjust(test_seqdenoise2,  censor_length = 5)
  test_seqdenoise2 = outseq(test_seqdenoise2, keep_flanks = TRUE)
  
  #no corrections applied for this one
  #test_seqdenoise2$outseq == t2
  expect_equal(test_seqdenoise2$outseq, t2)
  
  sweet_spot = "ACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCTTTAAGTTTATTAATTCGACTTGAATTAAGAACAATAAAAAATTTAATTGGTAATGATCAAATTTATAATGTAATTGTTACAGCCCACGCATTTATTATAATTTTTTTTATAGTTATACCAATTTTAATTGGAGGATTTGGAAATTGATTAGTCCCTATTATATTAGGAGCCCCAGATATAGCTTTTCCACGAATAAATAATATAAGATTCTGAATATTACCCCCTTCTTTATCTTTACTATCAATTAGAAGAATAGTAGAAACTGGGACCGGAACAGGATGAACTGTTTACCCACCTCTTTCTTCTGTAATTGCACATACTGGATCTTCAGTAGATTTTTCAATTTTTTCTCTACATATTGCAGGAATTTCATCTATTTTAGGAGCAATTAATTTTATTTCTACAATAATAAATATAAAAATTGAAAAATCTACTGAAG"
  test_seqdenoise2_trimmed = outseq(test_seqdenoise2, keep_flanks = FALSE)
  
  #test_seqdenoise2_trimmed$outseq == sweet_spot
  expect_equal(test_seqdenoise2_trimmed$outseq, sweet_spot)
  
 
  ############################################################################
  # Real read test 3
  ############################################################################
  # This one encounters a large run of ones and a large run of 0s about halfway through the barcode reg 
  
  t3 = 'AGGTCAACAAATCATAAAGATATTGGAACATTATATTTTTTATTTGGAATTTGAGCAGGAATAATTGGAACATCATTAAGATTATTAATTCGAATAGAATTAGGAAATCCTGGATCCTTAATTGGAGATGACCAAATTTATAATAACTATTGTAAACAGCCCATGCATTTATTATAATTTTTTTTATAGTAATACCTATTATAATTGGAGGATTTGGAAATTGATTAATTCCTTTAATATTAGGTGCACCAGATATAGCATTCCCCCGAATAAATAATATAAGATTTTGATTATTACCCCCTTCCATTTTTCTTTTAATTTCAAGAAGAATCGTAGAAAAATGGAGCAGGAACAGGATGAACTGTTTACCCCCCCTTATCTTCTAACACCGCTCATAGAGGAAGATCCGTAGATTTAGCCATTTTTTCTCTTCATTTAGCTGGAATTTCCTCAATTCTAGGAGCAGTAAATTTTATTTCTACAGTAATTAATATACGAGCTAAAAAAATAATATTTGACCAAATACCCCTATTTATTTGAGCTGTAGCTATTACTGCATTATTATTATTATTATCATTACCAGTTTTAGCAGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACATCTTTTTTTGACCCAGCTGGAGGAGGAGACCCAATTTTATACCAACATTTATTTTGATTTTTTGGACACCCTGAAGTTTA'

  #Below is behaving well despite initial doubt. The path encounters a large run of 0000s in the middle of the sequence and therefore terminates. Only possible improvement would be to add lots of Ns and pick up on the second run of 1s
  #possibly - could add a stop/reject condition if the adjust_seq length is less than the min_match?  
  test_seqdenoise3 = DNAseq(t3, name = "test3")
  test_seqdenoise3 = frame(test_seqdenoise3)
  # test_seqdenoise3$data$path
  # test_seqdenoise3$data$raw_removed_front
  # test_seqdenoise3$frame_dat
  
  test_seqdenoise3 = adjust(test_seqdenoise3, censor_length = 5)
  # test_seqdenoise3$adjusted_sequence
  test_seqdenoise3 = outseq(test_seqdenoise3)
  
  trimmed_and_framed3 = "NCATTATATTTTTTATTTGGAATTTGAGCAGGAATAATTGGAACATCATTAAGATTATTAATTCGAATAGAATTAGGAAATCCTGGATCCTTAATTGGAGATGACCAAATTTATAA"
  test_seqdenoise3_trimmed =  outseq(test_seqdenoise3, keep_flanks = FALSE)
  
  #test_seqdenoise3_trimmed$outseq == trimmed_and_framed3
  expect_equal(test_seqdenoise3_trimmed$outseq, trimmed_and_framed3)
  
  right_keep_expected =   "NCATTATATTTTTTATTTGGAATTTGAGCAGGAATAATTGGAACATCATTAAGATTATTAATTCGAATAGAATTAGGAAATCCTGGATCCTTAATTGGAGATGACCAAATTTATAATAACTATTGTAAACAGCCCATGCATTTATTATAATTTTTTTTATAGTAATACCTATTATAATTGGAGGATTTGGAAATTGATTAATTCCTTTAATATTAGGTGCACCAGATATAGCATTCCCCCGAATAAATAATATAAGATTTTGATTATTACCCCCTTCCATTTTTCTTTTAATTTCAAGAAGAATCGTAGAAAAATGGAGCAGGAACAGGATGAACTGTTTACCCCCCCTTATCTTCTAACACCGCTCATAGAGGAAGATCCGTAGATTTAGCCATTTTTTCTCTTCATTTAGCTGGAATTTCCTCAATTCTAGGAGCAGTAAATTTTATTTCTACAGTAATTAATATACGAGCTAAAAAAATAATATTTGACCAAATACCCCTATTTATTTGAGCTGTAGCTATTACTGCATTATTATTATTATTATCATTACCAGTTTTAGCAGGAGCTATTACTATATTATTAACAGATCGAAATTTAAATACATCTTTTTTTGACCCAGCTGGAGGAGGAGACCCAATTTTATACCAACATTTATTTTGATTTTTTGGACACCCTGAAGTTTA"

  test_seqdenoise3_keep_right = outseq(test_seqdenoise3, keep_flanks = 'right')

  expect_equal(test_seqdenoise3_keep_right$outseq, right_keep_expected)
  
  ############################################################################
  # Real read test 4
  ############################################################################
  
  #This one introduces an error into the output with a bunch of NAs
  #reverse compliment again the correct orientation
  t4 = 'GGTAGCTGTCAGAGTAGCTCGTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCAGGGTGTCCAAAAAATCAAAATAAATGTTGATATAAAATAGGGTCTCCCCCACCAGCTGGGTCAAAAAATGAGGTATTTAAATTTCGATCAGTTAATAATATAGTAATTGTTCCAGCTAAAACTGGTAAAGATAATAATAAAAGAAATGCAGTATACCATCTGCTCAAACAAATAATGGTATTTGATCAAATGATAAATTATTTAAACGTATATTAATAATTGTAGAAATAAAATTAATAGCTCCTAAAATAGATGAAATCCCAGCTAAATGAAGAGAAAAAATAGCTAAATCAACAGATCTTCCTCCATGAGCAATATTAGAAGAAAGAGGAGGATAAACTGTTCATCCAGTTCCTGCTCCATTTTCTACAATTCTACTTGAAATTAATAAAGTTAATGATGGGGGTAGAAGTCAAAAACTTATATTATTTATTCGGGGGAAGCTATATCTGGGGCTCCTAATATTAATGGTACTAATCAATTACCAAAT'
  
  test_seqdenoise4 = DNAseq(t4, name = "test4")
  test_seqdenoise4 = frame(test_seqdenoise4)
  # test_seqdenoise4$data$path
  # test_seqdenoise4$data$raw_removed_front
  # test_seqdenoise4$frame_dat
  
  test_seqdenoise4 = adjust(test_seqdenoise4,  censor_length = 5)
  # test_seqdenoise4$adjusted_sequence
  test_seqdenoise4 = outseq(test_seqdenoise4)
  
  #test_seqdenoise4$adjustment_count ==2
  expect_equal(test_seqdenoise4$adjustment_count, 2)
  

  #This one is an instance where the reverse compliment is the sequence that matches the path
  #therefore if we override the direction check, the sequence should be rejected
  test_seqdenoise4_trimmed = DNAseq(t4, name = "test4")
  test_seqdenoise4_trimmed = frame(test_seqdenoise4_trimmed, dir_check = FALSE)
  
  #test_seqdenoise4_trimmed$reject == TRUE
  expect_equal(test_seqdenoise4_trimmed$reject, TRUE)
  
  
  ############################################################################
  # Real read test 5
  ############################################################################
  
  #This one looks like its making legit corrections, double check and make sure the path and the output are reflective of one another
  t5 = 'ATTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGAATATTAGCCCCCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGAGCAGGAAACAGGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTCTTTAGCCAGTACTAGCTGGAGCTATTACTATAGCTACTTACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTATTTTGATTTTTTGGACATCCTGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACGTCTCTCGTCTGTGCCTACC'

  test_seqdenoise5 = DNAseq(t5, name = "test5")
  test_seqdenoise5 = frame(test_seqdenoise5)
  #path has 4 inserts scattered throughout a string of ones, plus a triple delete. lots of missing front and added end.
  #see if moving the back to the front reveals a circle?
  #test_seqdenoise5$data$path

  test_seqdenoise5 = adjust(test_seqdenoise5,  censor_length = 5)

  test_seqdenoise5 = outseq(test_seqdenoise5, adjust_limit  = 5)
  
  #based on path should apply 4 adjustments to the sequence
  #test_seqdenoise5$adjustment_count == 4
  expect_equal(test_seqdenoise5$adjustment_count, 4)
  
  #the output sequence, with edges reattached and the 
  #masks applied
  masked_output =  "ATTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGNNNNNNNNNNCCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGNNNNNNNNNNGGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTNNNNNNNNNNACTAGCTGGAGCTATTANNNNNNNNNNTACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTATTTTGATTTTTTGGACATCCTGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACGTCTCTCGTCTGTGCCTACC"
  #test_seqdenoise5$outseq == masked_output
  expect_equal(test_seqdenoise5$outseq, masked_output)
  
  #four inserts were detected and removed, so final length should be 4 less than the original
  #(nchar(test_seqdenoise5$outseq)+4) == nchar(t5)
  expect_equal((nchar(test_seqdenoise5$outseq)+4), nchar(t5))
  
  #should turn it into Ns, adjusted sequence is just the frame front
  #has the front edge added to the sequence.
  trimmed_and_framed5 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGNNNNNNNNNNCCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGNNNNNNNNNNGGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTNNNNNNNNNNACTAGCTGGAGCTATTANNNNNNNNNNTACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTA"
  test_seqdenoise5_trimmed =  outseq(test_seqdenoise5, keep_flanks = FALSE)
  
  #test_seqdenoise5_trimmed$outseq == trimmed_and_framed5
  expect_equal(test_seqdenoise5_trimmed$outseq, trimmed_and_framed5)
  
  #the trimmed output should be 654 in length, because of the triple delete in the path
  #nchar(test_seqdenoise5_trimmed$outseq) == 654
  expect_equal(nchar(test_seqdenoise5_trimmed$outseq), 654)
  
  
  ############################################################################
  # Real read test 6
  ############################################################################
  
  #this had NAs in the output?
  t6 = "GGTAGCGATCAGCTGAGCGCGTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCAGGATGACCAAAAAATCAAAATAAATGTTGATATAAAATAGGGTCTCCTCCTCCAGCAGGGTCAAAAAATTAAGTATTTAAATTTCGGTCTGTTAATAGTATAGTAATAGCTCCTGCTTAAACTGGTAAGGATAAAAGTAATAATAAAGCTGTAATAACAACTGATCAAACAAATAATGGTATTCGATCATATGAAATTCCATTAGATCGTATATTAATTACAGTAGTAATAAAATTTACTGCTCCTAAAATTCATGATTTTCATGAGATAAATATCACCTTATTCTCAATCCATTGTACCAACAAAAATTCTTTCCTATCTTCTATTCCTTCTAACTAATGCCGTGCAACTTGCATTTAAAAGACCTTAGGTCTTACAGGCGTATAGTTAAAATTATTATTATTATTATTATTAAAATAGATGATATACCTGCTTAATGAAGAGAAAAAATTGCTAAATCAACAGAGGCTCCTCTATGAGCAATTCTAGAAGAAAGAGGGGGATAAACTGTTCATCCTGTTCCAGCTCCAGTTTTCAACTATACTACTTATTAGTAATAAAGTTAAAGAAGGAGGTAATAGTCAAAACCTTATATTATTTATTCGAGGAAATGCTATATCTGGAGCTCCCAATATTAGAGGTACTAATCAATTTCCAAATCCTCCAATTATAATTGGTATTACCTATAAAAAAAATTATAACAAAAGCGTGAGCTGTAACAATTACATTATAAATTTGATCATCTCCAATCAAAGCTCCTGGATGACCAAGTTCTGCACGAATTAAAATTCTTAAAGAAGTTCCTACTATTCCAGCTCAAGTTCCATATAAAAAATATAATGTTCCAATATCTTTATGATTTGTTGAACTGTGACCTGAGTCAGCTACATGTTCGACTGCAGAGACTGCGACGAGACTACC"
  #this one the reverse compliment is the better match, and appears to be an ~420 run of good read with a single indel.
  #this one would also be solved by a minimum length on the adjust seq.
  
  test_seqdenoise6 = DNAseq(t6, name = "test6")
  test_seqdenoise6 = frame(test_seqdenoise6)
  #test_seqdenoise6$data$path
  
  test_seqdenoise6 = adjust(test_seqdenoise6,  censor_length = 5)
  
  test_seqdenoise6 = outseq(test_seqdenoise6, adjust_limit  = 5)
  
  #test_seqdenoise6$adjustment_count == 2
  expect_equal(test_seqdenoise6$adjustment_count, 2)
  
  
  ############################################################################
  # Real read test 7
  ############################################################################

  #This one has decent phred scores but appears to be getting completely censored.
  #^ is this one of the rev compliment reads? YES!
  t7 = 'GTGCGTGTGACTTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGGTGTCCAAAAAATCAAAATAAATGTTGATATAAAATTGGATCACCTCCTCCTGAAGGATCAAAAAATGAAGTATTTAAATTTCGATCAAATAATAATATTGTAATAGCTCCTGCTAATACAGGTAATGACAATAATAATAAAATAGCCGTTAATAATAATGATCAACAAAATAATCTTACTTTTTCTATTTTATATATTCTCATATTTAAAATAGTAGTAATAAAATTAATTGATGCTATAATTGAAGATATACCAGCAATATGCAATGAAAAGATTGATAAATCTACCGATGCTCCAGAATGAGATAAATTTATAGACAATGGTGGATATACAGTTCATCCTGTTCCAGAACCTATTCCAATAAATATTCTTGATATTAATAATATAATTCTTGGAGGTAATAATCAAAATCTTATATTATTTATTCGAGGAAATGCTATATCAGGAACCGATAATATTATAGGAACTAAAAAATTTCCAAAACCCCCTATTATTACTGGTATCACAAAAAAAAAAAATTATTACAAAAGCATGAATAGTTACAATAGCATTATAAATTTGATCATTTCCTAAAAATGAACCAGGATTTCCTAATTCTAAACGAATTATTATTCTTAATGATAATCCTATTATTCCTGCTCAAAATCCAAAAATAAAATATAAAATTCCAATATCTTTATGATTGGTTGAATGTGACCTGAGTCAGCTACATGTTCGACTGCTGCGTGAGAGAGAGACCTACCGTACATGTGCAGTCGAACATGTAGCTGACTCAGGTCACGTTCAACAAATCATAAAGTATAATGTTCCAATATCTTTATGATTTGTTGAACGTGACCTGAGTCAGCTACATGTTCGACTGCACATGTACATCAGCTCCTACC'

  test_seqdenoise7 = DNAseq(t7, name = "test7")
  test_seqdenoise7 = frame(test_seqdenoise7)

  test_seqdenoise7 = adjust(test_seqdenoise7)
  test_seqdenoise7$adjusted_sequence
  test_seqdenoise7 = outseq(test_seqdenoise7)
  
  #test_seqdenoise7$adjustment_count == 1
  expect_equal(test_seqdenoise7$adjustment_count, 1)
  
  })