
test_that("The denoising in a batch setting proceeds properly.", {

  #this string features a reported double insert to be covered up at 
  #positions 457 and 458 in the path.
  test_seq1 = 'GGTAGGAGCGTGTATACAGCGGCAGTCGAACATGTAGCTGACTCAGGTCACATTCAACAAATCATAAAGATATTGGAACATTATATTTTATTTTTGGAGCTTGATCTGGTATAGTAGGAACATCTTTAAGAATACTTATCCGAGCAGAATTAGGTCATCCAGGTACTTTTATTGGAGATGACCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTATCATAATTTTTTTTATAGTTATGCCTATCTTAATTGGAGGATTTGGAAATTGATTATTACCATTAATACTTGGGGCTCCAGATATAGCCTTCCCCCGAATAAATAATATAAGTTTTTGACTTCTTCCACCATCTTTAACTCTTCTTCTTTCTAGTTCCATTGTAGAAAATGGAGCAGGAACAGGATGAACTGTTTATCCTCCTCTTTCAGCAAGAATTGCACACAGAGGGGCTTCAGTAGATTTAGCTATCTTTTCTCTTCATTTAGCCGGAATCTCATCAATTTTAGGGTCAGTAAATTTTATTACAACAGGCAATTAATATACGAGCAAATGGGATTACTTTAGATCGAATACCTTTATTTGTATGATCAGTAGTAATTACAACAGTTCTTCTTCTTTTATCTCTTCCCGTATTAGCCGGAGCAATTACAATATTACTAACAGACCGAAATTTAAATACATCATTCTTTGACCCTGCAGGAGGAGGTGATCCAATTCTTTATCAACATTTATTTTGATTTTTTGGTACATCCATGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCAATTCTTTATCAACATTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCATGAGTGACGTGTAGCGCTACC'
  
  dat1 = DNAseq(test_seq1, name = 'testseq1' )
  
  dat1 = frame(dat1)
  
  names(dat1)
  
  dat1 = adjust(dat1, censor_length = 5)
  
  names(dat1)
  dat1$adjusted_sequence
  
  
  
  x = dat1
  censor_length = 5
  x$adjusted_sequence = adj_seq(x$frame_dat, x$data$path, censor_length = censor_length)
  
  frame_dat = x$frame_dat
  path_out = x$data$path
  censor_length = 5
  
  new_seq = c()
  #start at the first positon in the path
  org_seq_pos = 1
  path_pos = frame_dat$path_start
  path_end = frame_dat$path_end
  new_pos = 1
  
  censor_0s = c() #these were deleted from true seq  
  censor_2s = c() #these were inserted into true seq
  match_seen = FALSE
  
  triple_inserts = triple_ins(path_out, path_pos, path_end)
  
  x = path_out 
  path_start = path_pos
  path_end = path_end
  
  
  
  # this one
  test_seq2 = 'cgtatgtatgtcgcgcagtcgaacatgtagctgactcaggtcacggtcaacaaatcataaagatattggaactttatattttatttttggaatttgagcaggaatagtaggaacttctttaagtttattaattcgacttga
  attaagaacaataaaaaatttaattggtaatgatcaaatttataatgtaattgttacagcccacgcatttattataattttttttatagttataccaattttaattggaggatttggaaattgattagtccctattatattaggagccccagatatagcttttccacgaataaataatataagattctgaatattacccccttctttatctttactatcaattagaagaatagtagaaactgggaccggaacaggatgaactgtttacccacctctttcttctgtaattgcacatactggatcttcagtagatttttcaattttttctctacatattgcaggaatttcatctattttaggagcaattaattttatttctacaataataaatataaaaattgaaaaatctactgaagatccagtatgtgcaattacagaagaaagaggtgggtaaacagttcatcctgttccggtcccagtttctactattcttctaattaatagtaaagataaagaagggggtaatattcagaatcttatattatttattcgtggaaaagctatatctggggctcctaatataatagggactaatcaatttccaaatcctccaattaaaattggtataactataaaaaaaaattataataaatgcgtgggctgtaacaattacattataaatttgatcattaccaattaaattttttattgttcttaattcaagtcgaattaataaacttaaagaagttcctactattcctgctcaaattccaaaaataaaatataaagttccaatatctttatgattggttgaatgtgacctgagtcagctacatgttcgactgcgcgacatacatacg'
  
  dat2 = DNAseq(test_seq2, name = 'testseq2' )
  
  dat2 = frame(dat2)
  
  names(dat2)
  
  dat2 = adjust(dat2, censor_length = 5)
  
  dat2$adjusted_sequence
  dat2$adjustment_count
  
  names(dat1)
  dat1$adjusted_sequence
  
  
  
})