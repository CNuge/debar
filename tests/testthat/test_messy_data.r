test_that("Messy SEQUEL sequences are processed as anticipated", {
  
  ############################################################################
  # Real read test 1
  ############################################################################
  #this has a strong run of barcode data, there is a single reported error in the sequence (deleted character)
  # blast of the sequence indicates it is a missing t at the homopolymer tract. ~149-157
  #the path appears to get the position right :)
  
  t1 = "gctcctaaaattgaagaaactcctgctaaatgtaaagaaaaaattgctaaatctactgaagcaccaccatgagcaataatagaagataaaggtggataaacagttcatccagtaccagccccattttctactatacttctagccattaatagagtcaatgaaggcggtaataatcagaaactcatattatttattcgaggaaatgctatatccggagctcctaatattaatggaactaatcaatttccaaatcctccaattataattggtattactataaaaaaaattattacgaaagcatgggcagtgacaattacattataaatttgatcatcaccaattaatgctcctggatgtcctaattcagctcgaattaagattcttagagaagttcctactattccagctcatgctccaaatataaagtacaaagttccaatatctttatgattggttgactgtgacctgagtcagctacatgttcgactgcatagacagtgcgctgtctatgcagtcgaacatgtagctgactcaggtcacggtcaacaaatcataaagatattggaactttgtactttatatttggagcatgagctggaatagtaggaacttctctaagaatcttaattcgagctgaattaggacatccaggagcattaattggtgatgatcaaatttataatgtaattgtcactgcccatgctttcgtaataatttttttatagtaataccaattataattggaggatttggaaattgattagttccattaatattaggagctccagatatagcatttcctcgaataaataatatgagtctctgattattaccgccttcattgactctattaatggctagaagtatagtagaaaatggggctggtactggatgaactgtttatccacctttatcttctattattgctcatggtggtgcttcagtagatttagcaattttttctttacatttagcaggagtttcttcaattttaggagc"
  t1a = rev_comp(t1)
  
  test_seqdenoise = DNAseq(t1, name = "test1")
  test_seqdenoise = frame(test_seqdenoise)
  test_seqdenoise$data$path
  test_seqdenoise = frame(test_seqdenoise)
  test_seqdenoise = adjust(test_seqdenoise)
  test_seqdenoise = outseq(test_seqdenoise)
  
  #
  test1_corrected = "GCTCCTAAAATTGAAGAAACTCCTGCTAAATGTAAAGAAAAAATTGCTAAATCTACTGAAGCACCACCATGAGCAATAATAGAAGATAAAGGTGGATAAACAGTTCATCCAGTACCAGCCCCATTTTCTACTATACTTCTAGCCATTAATAGAGTCAATGAAGGCGGTAATAATCAGAAACTCATATTATTTATTCGAGGAAATGCTATATCCGGAGCTCCTAATATTAATGGAACTAATCAATTTCCAAATCCTCCAATTATAATTGGTATTACTATAAAAAAAATTATTACGAAAGCATGGGCAGTGACAATTACATTATAAATTTGATCATCACCAATTAATGCTCCTGGATGTCCTAATTCAGCTCGAATTAAGATTCTTAGAGAAGTTCCTACTATTCCAGCTCATGCTCCAAATATAAAGTACAAAGTTCCAATATCTTTATGATTGGTTGACTGTGACCTGAGTCAGCTACATGTTCGACTGCATAGACAGTGCGCTGTCTATGCAGTCGAACATGTAGCTGACTCAGGTCACGGTCAACAAATCATAAAGATATTGGAACTTTGTACTTTATATTTGGAGCATGAGCTGGAATAGTAGGAACTTCTCTAAGAATCTTAATTCGAGCTGAATTAGGACATCCAGGAGCATTAATTGGTGATGATCAAATTTATAATGTAATTGTCACTGCCCATGCTTTCGTAATAATTTNNNNNNNNNNNATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTTCCATTAATATTAGGAGCTCCAGATATAGCATTTCCTCGAATAAATAATATGAGTCTCTGATTATTACCGCCTTCATTGACTCTATTAATGGCTAGAAGTATAGTAGAAAATGGGGCTGGTACTGGATGAACTGTTTATCCACCTTTATCTTCTATTATTGCTCATGGTGGTGCTTCAGTAGATTTAGCAATTTTTTCTTTACATTTAGCAGGAGTTTCTTCAATTTTAGGAGC"

  #test_seqdenoise$outseq == test1_corrected
  expect_equal(test_seqdenoise$outseq, t1)
  
  
  ############################################################################
  # Real read test 2
  ############################################################################
  
  #issue here, the 296 characters outside of the denoised area are dropped
  #i.e. after the run of 2s happens, all of the following data is omitted in the subsequent output.
  #need adjust seq to keep these parts on the edges, find the logic that ceases the iteration if 
  # a large string of 1s or twos is encountered... consider how to proceed given that the later
  #example, t3 is benefitting from the trunction at this early point.
  t2 = "CGTATGTATGTCGCGCAGTCGAACATGTAGCTGACTCAGGTCACGGTCAACAAATCATAAAGATATTGGAACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCTTTAAGTTTATTAATTCGACTTGAATTAAGAACAATAAAAAATTTAATTGGTAATGATCAAATTTATAATGTAATTGTTACAGCCCACGCATTTATTATAATTTTTTTTATAGTTATACCAATTTTAATTGGAGGATTTGGAAATTGATTAGTCCCTATTATATTAGGAGCCCCAGATATAGCTTTTCCACGAATAAATAATATAAGATTCTGAATATTACCCCCTTCTTTATCTTTACTATCAATTAGAAGAATAGTAGAAACTGGGACCGGAACAGGATGAACTGTTTACCCACCTCTTTCTTCTGTAATTGCACATACTGGATCTTCAGTAGATTTTTCAATTTTTTCTCTACATATTGCAGGAATTTCATCTATTTTAGGAGCAATTAATTTTATTTCTACAATAATAAATATAAAAATTGAAAAATCTACTGAAGATCCAGTATGTGCAATTACAGAAGAAAGAGGTGGGTAAACAGTTCATCCTGTTCCGGTCCCAGTTTCTACTATTCTTCTAATTAATAGTAAAGATAAAGAAGGGGGTAATATTCAGAATCTTATATTATTTATTCGTGGAAAAGCTATATCTGGGGCTCCTAATATAATAGGGACTAATCAATTTCCAAATCCTCCAATTAAAATTGGTATAACTATAAAAAAAAATTATAATAAATGCGTGGGCTGTAACAATTACATTATAAATTTGATCATTACCAATTAAATTTTTTATTGTTCTTAATTCAAGTCGAATTAATAAACTTAAAGAAGTTCCTACTATTCCTGCTCAAATTCCAAAAATAAAATATAAAGTTCCAATATCTTTATGATTGGTTGAATGTGACCTGAGTCAGCTACATGTTCGACTGCGCGACATACATACG"
  t2a = rev_comp(t2)
  
  test_seqdenoise2 = DNAseq(t2, name = "test1")
  test_seqdenoise2 = frame(test_seqdenoise2)
  test_seqdenoise2 = adjust(test_seqdenoise2)
  test_seqdenoise2 = outseq(test_seqdenoise2)
  
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
  
  test_seqdenoise3 = DNAseq(t3, name = "test3")
  test_seqdenoise3 = frame(test_seqdenoise3)
  # test_seqdenoise3$data$path
  # test_seqdenoise3$data$raw_removed_front
  # test_seqdenoise3$frame_dat
  
  test_seqdenoise3 = adjust(test_seqdenoise3)
  # test_seqdenoise3$adjusted_sequence
  test_seqdenoise3 = outseq(test_seqdenoise3)
  
  trimmed_and_framed3 = "NCATTATATTTTTTATTTGGAATTTGAGCAGGAATAATTGGAACATCATTAAGATTATTAATTCGAATAGAATTAGGAAATCCTGGATCCTTAATTGGAGATGACCAAATTTATAA"
  test_seqdenoise3_trimmed =  outseq(test_seqdenoise3, keep_flanks = FALSE)
  
  #test_seqdenoise3_trimmed$outseq == trimmed_and_framed3
  expect_equal(test_seqdenoise3_trimmed$outseq, trimmed_and_framed3)
  
  ############################################################################
  # Real read test 4
  ############################################################################
  
  #This one introduces an error into the output with a bunch of NAs
  t4 = 'GGTAGCTGTCAGAGTAGCTCGTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCAGGGTGTCCAAAAAATCAAAATAAATGTTGATATAAAATAGGGTCTCCCCCACCAGCTGGGTCAAAAAATGAGGTATTTAAATTTCGATCAGTTAATAATATAGTAATTGTTCCAGCTAAAACTGGTAAAGATAATAATAAAAGAAATGCAGTATACCATCTGCTCAAACAAATAATGGTATTTGATCAAATGATAAATTATTTAAACGTATATTAATAATTGTAGAAATAAAATTAATAGCTCCTAAAATAGATGAAATCCCAGCTAAATGAAGAGAAAAAATAGCTAAATCAACAGATCTTCCTCCATGAGCAATATTAGAAGAAAGAGGAGGATAAACTGTTCATCCAGTTCCTGCTCCATTTTCTACAATTCTACTTGAAATTAATAAAGTTAATGATGGGGGTAGAAGTCAAAAACTTATATTATTTATTCGGGGGAAGCTATATCTGGGGCTCCTAATATTAATGGTACTAATCAATTACCAAAT'
  
  test_seqdenoise4 = DNAseq(t4, name = "test4")
  test_seqdenoise4 = frame(test_seqdenoise4)
  # test_seqdenoise4$data$path
  # test_seqdenoise4$data$raw_removed_front
  # test_seqdenoise4$frame_dat
  
  test_seqdenoise4 = adjust(test_seqdenoise4)
  # test_seqdenoise4$adjusted_sequence
  test_seqdenoise4 = outseq(test_seqdenoise4)
  
  #should turn it into Ns, adjusted sequence is just the frame front
  trimmed_and_framed4 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
  test_seqdenoise4_trimmed =  outseq(test_seqdenoise4, keep_flanks = FALSE)
 
  #test_seqdenoise4_trimmed$outseq == trimmed_and_framed4
  expect_equal(test_seqdenoise4_trimmed$outseq, trimmed_and_framed4)
  
  
  
  ############################################################################
  # Real read test 5
  ############################################################################
  
  #This one looks like its making legit corrections, double check and make sure the path and the output are reflective of one another
  t5 = 'ATTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGAATATTAGCCCCCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGAGCAGGAAACAGGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTCTTTAGCCAGTACTAGCTGGAGCTATTACTATAGCTACTTACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTATTTTGATTTTTTGGACATCCTGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACGTCTCTCGTCTGTGCCTACC'

  test_seqdenoise5 = DNAseq(t5, name = "test4")
  test_seqdenoise5 = frame(test_seqdenoise5)
  #path has 3 inserts scattered throughout a string of ones, plus a triple delete. lots of missing front and added end.
  #see if moving the back to the front reveals a circle?
  #test_seqdenoise5$data$path

  test_seqdenoise5 = adjust(test_seqdenoise5)

  test_seqdenoise5 = outseq(test_seqdenoise5, adjust_limit  = 5)
  
  #based on path should apply 4 adjustments to the sequence
  #test_seqdenoise5$adjustment_count == 4
  expect_equal(test_seqdenoise5$adjustment_count, 4)
  
  #the output sequence, with edges reattached and the 
  #masks applied
  masked_output =  "ATTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGANNNNNNNNNNCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGANNNNNNNNNNGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTCNNNNNNNNNNCTAGCTGGAGCTATTACNNNNNNNNNNACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTATTTTGATTTTTTGGACATCCTGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACGTCTCTCGTCTGTGCCTACC"
  #test_seqdenoise5$outseq == masked_output
  expect_equal(test_seqdenoise5$outseq, masked_output)
  
  #four inserts were detected and removed, so final length should be 4 less than the original
  #(nchar(test_seqdenoise5$outseq)+4) == nchar(t5)
  expect_equal((nchar(test_seqdenoise5$outseq)+4), nchar(t5))
  
  #should turn it into Ns, adjusted sequence is just the frame front
  #has the front edge added to the sequence.
  trimmed_and_framed4 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTGACTAGTTCCTTTAAAATTAGGAGCCCCAGATATAGCATTTCCTCGGATAAATAATATAAGTTTTTGANNNNNNNNNNCTTCATTAACTTTACTTTTATCAAGCTCTATTGTAGAAAATGGANNNNNNNNNNGTTGAACTGTTTACCCTCCTTTATCTTCTGGGATTGCCCATGCAGGAGCTTCTGTTGATTTAGCTATTTTTTCTCTTCACTTAGCTGAAATTTCCTCAATTTTAGGAGCAGTAAATTTTATTACAACTGTAATTAATATACGGTCTAGAGGAATTACTTTAGATCGAATACCTTTATTTGTTTGATCAGTAGTAATTACTGCAATTTTATTACTTCTCTCNNNNNNNNNNCTAGCTGGAGCTATTACNNNNNNNNNNACAGACCGAAATTTAAATACCTCTTTTTTTGACCCAGCAGGGGGAGGGGACCCAATTCTTTATCAACACTTA"
  test_seqdenoise5_trimmed =  outseq(test_seqdenoise5, keep_flanks = FALSE)
  
  #test_seqdenoise5_trimmed$outseq == trimmed_and_framed4
  expect_equal(test_seqdenoise5_trimmed$outseq, trimmed_and_framed4)
  
  #the trimmed output should be 654 in length, because of the triple delete in the path
  #nchar(test_seqdenoise5_trimmed$outseq)== 654
  expect_equal(nchar(test_seqdenoise5_trimmed$outseq), 654)
  
  
  
  ############################################################################
  # Real read test 6
  ############################################################################
  
  #this had NAs in the output?
  t6 = "GGTAGCGATCAGCTGAGCGCGTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCAGGATGACCAAAAAATCAAAATAAATGTTGATATAAAATAGGGTCTCCTCCTCCAGCAGGGTCAAAAAATTAAGTATTTAAATTTCGGTCTGTTAATAGTATAGTAATAGCTCCTGCTTAAACTGGTAAGGATAAAAGTAATAATAAAGCTGTAATAACAACTGATCAAACAAATAATGGTATTCGATCATATGAAATTCCATTAGATCGTATATTAATTACAGTAGTAATAAAATTTACTGCTCCTAAAATTCATGATTTTCATGAGATAAATATCACCTTATTCTCAATCCATTGTACCAACAAAAATTCTTTCCTATCTTCTATTCCTTCTAACTAATGCCGTGCAACTTGCATTTAAAAGACCTTAGGTCTTACAGGCGTATAGTTAAAATTATTATTATTATTATTATTAAAATAGATGATATACCTGCTTAATGAAGAGAAAAAATTGCTAAATCAACAGAGGCTCCTCTATGAGCAATTCTAGAAGAAAGAGGGGGATAAACTGTTCATCCTGTTCCAGCTCCAGTTTTCAACTATACTACTTATTAGTAATAAAGTTAAAGAAGGAGGTAATAGTCAAAACCTTATATTATTTATTCGAGGAAATGCTATATCTGGAGCTCCCAATATTAGAGGTACTAATCAATTTCCAAATCCTCCAATTATAATTGGTATTACCTATAAAAAAAATTATAACAAAAGCGTGAGCTGTAACAATTACATTATAAATTTGATCATCTCCAATCAAAGCTCCTGGATGACCAAGTTCTGCACGAATTAAAATTCTTAAAGAAGTTCCTACTATTCCAGCTCAAGTTCCATATAAAAAATATAATGTTCCAATATCTTTATGATTTGTTGAACTGTGACCTGAGTCAGCTACATGTTCGACTGCAGAGACTGCGACGAGACTACC"
  t6a = rev_comp(t6)
  #this one the reverse compliment is the better match, and appears to be an ~420 run of good read with a single indel.
  #
  
  test_seqdenoise6 = DNAseq(t6a, name = "test6")
  test_seqdenoise6 = frame(test_seqdenoise6)
  #test_seqdenoise6$data$path
  
  test_seqdenoise6 = adjust(test_seqdenoise6)
  
  test_seqdenoise6 = outseq(test_seqdenoise6, adjust_limit  = 5)
  
  "GGTAGCGATCAGCTGAGCGCGTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCAGGATGACCAAAAAATCAAAATAAATGTTGATATAAAATAGGGTCTCCTCCTCCAGCAGGGTCAAAAAATTAAGTATTTAAATTTCGGTCTGTTAATAGTATAGTAATAGCTCCTGCTTAAACTGGTAAGGATAAAAGTAATAATAAAGCTGTAATAACAACTGATCAAACAAATAATGGTATTCGATCANANANANANANANANANANANANANANANANANANANANANANANANANANANANANANANANANNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTATGAAATTCCATTAGATCGTATATTAATTACAGTAGTAATAAAATTTACTGCTCCTAAAATTCATGATTTTCATGAGATAAATATCACCTTATTCTCAATCCATTGTACCAACAAAAATTCTTTCCTATCTTCTATTCCTTCTAACTAATGCCGTGCAACTTGCATTTAAAAGACCTTAGGTCTTACAGGCGTATAGTTAAAATTATTATTATTATTATTATTAAAATAGATGATATACCTGCTTAATGAAGAGAAAAAATTGCTAAATCAACAGAGGCTCCTCTATGAGCAATTCTAGAAGAAAGAGGGGGATAAACTGTTCATCCTGTTCCAGCTCCAGTTTTCAACTATACTACTTATTAGTAATAAAGTTAAAGAAGGAGGTAATAGTCAAAACCTTATATTATTTATTCGAGGAAATGCTATATCTGGAGCTCCCAATATTAGAGGTACTAATCAATTTCCAAATCCTCCAATTATAATTGGTATTACCTATAAAAAAAATTATAACAAAAGCGTGAGCTGTAACAATTACATTATAAATTTGATCATCTCCAATCAAAGCTCCTGGATGACCAAGTTCTGCACGAATTAAAATTCTTAAAGAAGTTCCTACTATTCCAGCTCAAGTTCCATATAAAAAATATAATGTTCCAATATCTTTATGATTTGTTGAACTGTGACCTGAGTCAGCTACATGTTCGACTGCAGAGACTGCGACGAGACTACC"
  x = test_seqdenoise6
  outseq = paste(c(x$data$raw_removed_front, #part removed in front of first frame() call
                   x$frame_dat$removed_lead, #part removed in second frame() call
                   x$adjusted_sequence[dashes_rm:length(x$adjusted_sequence)], #the 'sweet spot' that is denoised
                   x$data$adjusted_trimmed, # part removed in adj_seq
                   x$frame_dat$removed_end, #part removed in second frame() call
                   x$data$raw_removed_end), # part removed from the end of first frame() call
                 collapse = "")
  
  ############################################################################
  # Real read test 7
  ############################################################################

  
  library(seqinr)
  t7a = c2s(rev(comp(s2c(t7))))
  
  #This one has decent phred scores but appears to be getting completely censored.
  #^ is this one of the rev compliment reads? YES!
  t7 = 'GTGCGTGTGACTTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGGTGTCCAAAAAATCAAAATAAATGTTGATATAAAATTGGATCACCTCCTCCTGAAGGATCAAAAAATGAAGTATTTAAATTTCGATCAAATAATAATATTGTAATAGCTCCTGCTAATACAGGTAATGACAATAATAATAAAATAGCCGTTAATAATAATGATCAACAAAATAATCTTACTTTTTCTATTTTATATATTCTCATATTTAAAATAGTAGTAATAAAATTAATTGATGCTATAATTGAAGATATACCAGCAATATGCAATGAAAAGATTGATAAATCTACCGATGCTCCAGAATGAGATAAATTTATAGACAATGGTGGATATACAGTTCATCCTGTTCCAGAACCTATTCCAATAAATATTCTTGATATTAATAATATAATTCTTGGAGGTAATAATCAAAATCTTATATTATTTATTCGAGGAAATGCTATATCAGGAACCGATAATATTATAGGAACTAAAAAATTTCCAAAACCCCCTATTATTACTGGTATCACAAAAAAAAAAAATTATTACAAAAGCATGAATAGTTACAATAGCATTATAAATTTGATCATTTCCTAAAAATGAACCAGGATTTCCTAATTCTAAACGAATTATTATTCTTAATGATAATCCTATTATTCCTGCTCAAAATCCAAAAATAAAATATAAAATTCCAATATCTTTATGATTGGTTGAATGTGACCTGAGTCAGCTACATGTTCGACTGCTGCGTGAGAGAGAGACCTACCGTACATGTGCAGTCGAACATGTAGCTGACTCAGGTCACGTTCAACAAATCATAAAGTATAATGTTCCAATATCTTTATGATTTGTTGAACGTGACCTGAGTCAGCTACATGTTCGACTGCACATGTACATCAGCTCCTACC'

  test_seqdenoise7 = DNAseq(t7, name = "test7")
  test_seqdenoise7 = frame(test_seqdenoise7)
  test_seqdenoise7$data$path
  test_seqdenoise7$data$raw_removed_front
  test_seqdenoise7$frame_dat
  
  test_seqdenoise7 = adjust(test_seqdenoise7)
  test_seqdenoise7$adjusted_sequence
  test_seqdenoise7 = outseq(test_seqdenoise7)
  
  #test_seqdenoise7$outseq == t7
  expect_equal(test_seqdenoise7$outseq, t7)
  
  
  #truncates to essentially nothing
  trimmed_and_framed7 = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGTGCGTGTGACTTG" 
  test_seqdenoise7_trimmed =  outseq(test_seqdenoise7, keep_flanks = FALSE)
  
  #test_seqdenoise7_trimmed$outseq == trimmed_and_framed7
  expect_equal(test_seqdenoise7_trimmed$outseq, trimmed_and_framed7)
  
  
  
  
  test_seqdenoise7a = DNAseq(t7a, name = "test7a")
  test_seqdenoise7a = frame(test_seqdenoise7a)
  test_seqdenoise7a$data$path
  test_seqdenoise7a$data$raw_removed_front
  test_seqdenoise7a$frame_dat

  test_seqdenoise7a = adjust(test_seqdenoise7a)
  test_seqdenoise7a$adjusted_sequence
  test_seqdenoise7a = outseq(test_seqdenoise7a)

  trimmed_and_framed7 =
  test_seqdenoise7a_trimmed =  outseq(test_seqdenoise7a, keep_flanks = FALSE)

  #test_seqdenoise7a_trimmed$outseq == trimmed_and_framed4
  expect_equal(test_seqdenoise7a_trimmed$outseq, trimmed_and_framed7)

})