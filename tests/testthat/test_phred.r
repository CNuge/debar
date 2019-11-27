test_that("fastq phred scores can be processed and maintained correctly", {


  ##############################################################3
  # specific examples
  # test 1 - no changes to the phred
  t1_seq = 'GGTAGGTCTTCTCTCACGCAGCAGTCGAACATGTAGCTGACTCAGGTCACGTTCAACAAATCATAAAGATATTGGTACATTATATTTTCTTTTTGGAATTTGAGCAGGAATAGTAGGAACATCACTAAGATTATTAATTCGTATAGAATTAAGAACAATTAGAAATTTAATTGGAAATGATCAAATTTATAATGTAATTGTAACTGCTCATGCTTTTATTATAATTTTCTTCATAGTAATACCAATTTTAATTGGAGGATTTGGAAATTGATTAATTCCAATTATACTAGGAGCCCCAGATATAGCTTTTCCTCGAATAAATAATATAAGATTTTGAATATTACCCCCATCTTTATCTTTATTATTAATTAGAAGAATAGTAGAAACTGGAACAGGAACAGGATGAACAGTTTACCCACCCTTATCATCAGTAATTGCTCATACAGGTTCATCTGTAGATTTCTCTATTTTTTCATTACATATTGCAGGTATTTCATCTATTTTAGGAGCTATTAATTTTATTTCAACAATAATAAATATAAAAATTAAATT'
  t1_phred = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      
  test_seqdenoise1 = DNAseq(t1_seq, name = "test1", phred = t1_phred)
  test_seqdenoise1 = frame(test_seqdenoise1)
  test_seqdenoise1 = adjust(test_seqdenoise1)
  test_seqdenoise1 = outseq(test_seqdenoise1)
  
  #test_seqdenoise1$phred== t1_phred
  expect_equal(test_seqdenoise1$phred, t1_phred)
  
  t2_seq = 'GGTAGAGATGTAGCACATCATTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGATGACCAAAAAATCAAAATAAATGTTGATATAAGACTGGATCTCCCCCTCCAGATGGATCAAAAAATGATGTATTTAAATTTCGATCAAATAATAATATAGTAAATAGCTCCAGCTAATACTGGTAATGATAATAATAATAAAATTGAAGTTAATAATATCGATCATCTAAATAATGTAACTTTATCTATTTTTGTAAATCTTATATTAATAATTGTTCTAATAAAAATTAATTGATGCTAAAATAGAAGAAATTCCTGCAATATGTAAAGAAAAAATAGATAAATCTACTGATATTCTTCTATGTGATATACTTAATGACAATAGTGGATATACAGTTCATCCAGTTCCAGTTCCAGTTCCAATAAATATTCTAAATATTAATATTATTATTCTAGGAATTAATAATCAAAATCTTATATTATTTATTCGTGGAAATGATATATCAGCAACATTTAATATTAAAGGAA'
  t2_phred = '~~~~~~~~~~~~z~~~~~~~R~~~~~~~~}~~~~~~~~{~~Hv~y~~}~~~}w~~~f~~~i~~~~~q~8~~~~~~~N~~~~e~~~~x~~~~~np~v~~q~~~~~~M~~~~~t~~~~~~~~~~B~~~~~s~~~~~~>~~y~~T~~~~~~~_~~~o~~~~~~~~~~~T~~~~~~~~~~~~o~~~~~v~~~~~~~~~~~~~~{~~U~~~~~~~~~{~~~~~~~~~~~~~~r~~~v~~~~~~~~b~~}~~~~~~~B~~~~~~t~~{~`~~~~~~~~~~~~~~`~~~s~~~2~~~~~~~v~~~~~~~R~~~~~~~~~d~~~~r~~~~~~~~~~~~~~~D~~~~~~~~~~b~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~z~~~~~~~~~L~~~~~~~~~y~~~~~l~~~u~~p~~~~p~~~I~~~}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~o~~~~~~~z~~~v~~~~~~o~~~~~~~~~~~~~~~~~~~~~~~~~~~k~~{~`~'
    
  test_seqdenoise2 = DNAseq(t2_seq, name = "test2", phred = t2_phred)
  test_seqdenoise2 = frame(test_seqdenoise2)
  test_seqdenoise2 = adjust(test_seqdenoise2)
  test_seqdenoise2 = outseq(test_seqdenoise2)
  
  #nchar(test_seqdenoise2$outphred) ==  (nchar(t2_phred) - 1)
  expect_equal(nchar(test_seqdenoise2$outphred), (nchar(t2_phred) - 1))
  
  #                                                        #after del here
  t3_seq = 'CACGAATAAATAATATAAGATTTTGATTTTTGCCCCCCTCATTAACTTACTGTTAATAAGAAGAGTTGTAGAAACAGGAACTGGGACAGGTTGAACAGTTTATCCTCCTTTATCATCAATTATTGCTCATACTGGAGCCTCTGTAGATTTATCTATTTTCTCTTTTCATATAGCCGGAATTTCATCCATTCTAGGAGCTATTAACTTTATTACAACTATAATTAATATACGAGTAAAAAATATTAAATTTGACCAAATCCCTTTATTCTCCTGATCAGTAATTATTACAGCTATTTTACTTTTATTATCTCTACCAGTCTTAGCCGGAGCTATTACTATATTATTAACAGATCGAAATCTTAATACCTCTTTTTTCGACCCCAGGGGAGGAGGGGATCCTATTTTATACCAACACTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACTCTGCTCTGACTCTCCTACC'
  t3_seq_out= "CACGAATAAATAATATAAGATTTTGATTTTTGCCCCCCTCATTAACTTNACTGTTAATAAGAAGAGTTGTAGAAACAGGAACTGGGACAGGTTGAACAGTTTATCCTCCTTTATCATCAATTATTGCTCATACTGGAGCCTCTGTAGATTTATCTATTTTCTCTTTTCATATAGCCGGAATTTCATCCATTCTAGGAGCTATTAACTTTATTACAACTATAATTAATATACGAGTAAAAAATATTAAATTTGACCAAATCCCTTTATTCTCCTGATCAGTAATTATTACAGCTATTTTACTTTTATTATCTCTACCAGTCTTAGCCGGAGCTATTACTATATTATTAACAGATCGAAATCTTAATACCTCTTTTTTCGACCCCAGGGGAGGAGGGGATCCTATTTTATACCAACACTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACTCTGCTCTGACTCTCCTACC"
  
  t3_phred = '~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
  
  t3_phred_out_expected = '~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~*~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~' 
                         
  test_seqdenoise3 = DNAseq(t3_seq, name = "test3", phred = t3_phred)
  
  #phred check should not flag this one
  test_seqdenoise3 = phred_check(test_seqdenoise3)
  #test_seqdenoise3$reject == FALSE
  expect_equal(test_seqdenoise3$reject, FALSE)
  
  test_seqdenoise3 = frame(test_seqdenoise3)
  test_seqdenoise3 = adjust(test_seqdenoise3, censor_length = 0)
  test_seqdenoise3 = outseq(test_seqdenoise3)
  
  #test_seqdenoise3$outseq == t3_seq_out
  expect_equal(test_seqdenoise3$outseq, t3_seq_out)
  #test_seqdenoise3$outphred == t3_phred_out_expected
  expect_equal(test_seqdenoise3$outphred, t3_phred_out_expected)
  
  
  ########################
  # This one has doctored low QV and should therefore be rejected
  t4_seq = "AGTTAACAAATCATAAAGATATTGGTACCTTATACTTTATTTTTGGAGCTTGAGCTGGTATAATTGGCACATCCCTTAGTATTATTATTCGTGCAGAATTAGGACATCCAGGTACATTTATTGGTAATGACNNNNNNNNNNNTGTTATTGTTACTGCTCATGCCTTTATTATAATTTTTTTCATAGTTATACCAATTATAATCGGTGGATTTGGAAATTGACTTGTACCTTTAATATTAGGTGCCCCCGATATAGCCTTTCCCCNNNNNNNNNNTATAAGATTTTGAATATTACCCCCTTCTATTACTTTATTACTAACAAGAAGACTTGTAGAGAATGGAGCAGGTACCGGATGAACAGTTTACCCACCATTAGCCGCTAATATTTCTCATGCTGGATCCTCAGTTGACCTAGCAATTTTCTCTTTACACTTAGCTGGGATTTCTTCCATTCTAGGAGCTGTAAATTTTATTACTACAATTATCAATATACGATCAAATGGAATTTCTTTCGACCGTATGCCTTTATTTGTATGATCTGTCTTTATTACAGCTATTTTATTACTTTTATCATTGCCTGTTTTAGCCGGAGCAATTACAATATTATTAACTGATCGTAATTTAAATACATCATTCTTTGACCCTGCAGGAGGGGGAGATCCAATTCTTTACCAACATTTATTTTGATTTTTGGTCATCCATGAAGTTTA"
  t4_phred="euq_9f[cx+>w[ng$qcyerornojy^:ebaqgvigh}q$nsaie]wkmlb7vmo=xfh{{4ezz@g$2>_>gpduqv$d8aslz=_kreydpogohj}shpfbx;gh_=i=l^pzu{|>[jppudcrmozuyypviked_h>oe5$0nw^qbd_&$axh_jdt0[us?btbdqbg[c9yzxeqrbtpk$^nysxcjau]c{zyt9jijefrjthkzsx[vnce;cqi<n-w>$]*x<lnl]q99jopiqhjas[qq2bvdpitly/gp]d|tni:winz{:ci}rw[hm[ub_ktmfw_evvdxwlouipinqjn_lisfsb)jcvrf]dhl_^$:w[]gfjb$wx?als{{1cj$kuxlhszp_qdr[:gex[la7$d_a$nbkxws0$.$5kvb;-plhaclh>^************************************************************************************************************************************************************************************************************************************************************************************************************"
  
  test_seqdenoise4 = DNAseq(t4_seq, name = "test4", phred = t4_phred)
  test_seqdenoise4 = phred_check(test_seqdenoise4)
  
  #test_seqdenoise4$reject == TRUE
  expect_equal(test_seqdenoise4$reject, TRUE)
  
  
  #####################
  # This one appears to be dropping the phred scores in the fq output, see if you can pin down the bug
  #@m54146_190221_235805/15859816/ccs
  t5_seq = "GAGTCAGAGATCTATGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGGTGTCCAAAAAAATCAAAATAAATGTTGATATAAAATTGGATCACCTCCTCCTGCTGGGTCAAAGAAAGATGTATTTAAATTTCGGTCGGTTAATAATATTGTAATAGCTCCAGCTAATACAGGTAGAGATAAAAGTAATAAAATAGCAGTAATTACAACAGATCATACAAATAAAGGTATTCGATCAAATGTAATACCTGTTGATCGTATATTAATTACAGTAGTAATGAAATTTACGGCTCCTAAAATAGAAGAAATACCAGCTAAATGTAATGAAAAAAATAGCTAAATCAACTGAGGCTCCTCCATGAGCAATCCCTGCTGATAGAGGAGGGTAAACTGTTCAACCAGTTCCAGCTCCATTTTCTACTATACTACTAGCTAATAATAATGTAAGAGAAGGGGGTAATATTCAAAAACTTATATTACCCCCTTCTCTTACATTATTATTAGCTAGTAGTATAGTAGAAAATGGAGCTGGAACTGGTTGAACAGTTTACCCTCCTCTATCAGCAGGGATTGCTCATGGAGGAGCCTCAGTTGATTTAGCTATTTTTTCATTACATTTAGCTGGTACTTCTTCTATTTTAGGAGCCGTAAATTTCATTACTACTGTAATTAATATACGATCAACAGGTATTACATTTGATCGAATACCTTTATTTGTATGATCTGTTGTAATTACTGCTATTTTATTACTTTTATCTCTACCTGTATTAGCTGGAGCTATTACAATATTATTAACCGACCGAAATTTAAATACATCTTTCTTTGACCCAGCAGGAGGAGGTGATCCAATTTTATATCAACATTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCATAGATCTCTGACTC"
  t5_phred = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~V~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~0~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~t~~~~~~~~~~~~a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~V~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  test_seqdenoise5 = DNAseq(t5_seq, name = "test5", phred = t5_phred)
  
  test_seqdenoise5 = frame(test_seqdenoise5)
  test_seqdenoise5 = adjust(test_seqdenoise5)
  test_seqdenoise5 = outseq(test_seqdenoise5)
  
  test_seqdenoise5$outseq
  test_seqdenoise5$outphred
})