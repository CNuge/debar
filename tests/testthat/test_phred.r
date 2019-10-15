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
  
  
  ##############################################################3
  # specific examples
  #
  t1_seq = 'GGTAGGTCTTCTCTCACGCAGCAGTCGAACATGTAGCTGACTCAGGTCACGTTCAACAAATCATAAAGATATTGGTACATTATATTTTCTTTTTGGAATTTGAGCAGGAATAGTAGGAACATCACTAAGATTATTAATTCGTATAGAATTAAGAACAATTAGAAATTTAATTGGAAATGATCAAATTTATAATGTAATTGTAACTGCTCATGCTTTTATTATAATTTTCTTCATAGTAATACCAATTTTAATTGGAGGATTTGGAAATTGATTAATTCCAATTATACTAGGAGCCCCAGATATAGCTTTTCCTCGAATAAATAATATAAGATTTTGAATATTACCCCCATCTTTATCTTTATTATTAATTAGAAGAATAGTAGAAACTGGAACAGGAACAGGATGAACAGTTTACCCACCCTTATCATCAGTAATTGCTCATACAGGTTCATCTGTAGATTTCTCTATTTTTTCATTACATATTGCAGGTATTTCATCTATTTTAGGAGCTATTAATTTTATTTCAACAATAATAAATATAAAAATTAAATT'
  t1_phred = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
      
  test_seqdenoise1 = DNAseq(t1_seq, name = "test1", phred = t1_phred)
  test_seqdenoise1 = frame(test_seqdenoise1)
  test_seqdenoise1 = adjust(test_seqdenoise1)
  test_seqdenoise1 = outseq(test_seqdenoise1)
  
  
  t2_seq = 'GGTAGAGATGTAGCACATCATTGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGATGACCAAAAAATCAAAATAAATGTTGATATAAGACTGGATCTCCCCCTCCAGATGGATCAAAAAATGATGTATTTAAATTTCGATCAAATAATAATATAGTAAATAGCTCCAGCTAATACTGGTAATGATAATAATAATAAAATTGAAGTTAATAATATCGATCATCTAAATAATGTAACTTTATCTATTTTTGTAAATCTTATATTAATAATTGTTCTAATAAAAATTAATTGATGCTAAAATAGAAGAAATTCCTGCAATATGTAAAGAAAAAATAGATAAATCTACTGATATTCTTCTATGTGATATACTTAATGACAATAGTGGATATACAGTTCATCCAGTTCCAGTTCCAGTTCCAATAAATATTCTAAATATTAATATTATTATTCTAGGAATTAATAATCAAAATCTTATATTATTTATTCGTGGAAATGATATATCAGCAACATTTAATATTAAAGGAA'
  t2_phred = '~~~~~~~~~~~~z~~~~~~~R~~~~~~~~}~~~~~~~~{~~Hv~y~~}~~~}w~~~f~~~i~~~~~q~8~~~~~~~N~~~~e~~~~x~~~~~np~v~~q~~~~~~M~~~~~t~~~~~~~~~~B~~~~~s~~~~~~>~~y~~T~~~~~~~_~~~o~~~~~~~~~~~T~~~~~~~~~~~~o~~~~~v~~~~~~~~~~~~~~{~~U~~~~~~~~~{~~~~~~~~~~~~~~r~~~v~~~~~~~~b~~}~~~~~~~B~~~~~~t~~{~`~~~~~~~~~~~~~~`~~~s~~~2~~~~~~~v~~~~~~~R~~~~~~~~~d~~~~r~~~~~~~~~~~~~~~D~~~~~~~~~~b~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~z~~~~~~~~~L~~~~~~~~~y~~~~~l~~~u~~p~~~~p~~~I~~~}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(~~~~~o~~~~~~~z~~~v~~~~~~o~~~~~~~~~~~~~~~~~~~~~~~~~~~k~~{~`~'
    
  test_seqdenoise2 = DNAseq(t2_seq, name = "test2", phred = t2_phred)
  test_seqdenoise2 = frame(test_seqdenoise2)
  test_seqdenoise2 = adjust(test_seqdenoise2)
  test_seqdenoise2 = outseq(test_seqdenoise2)
  
  t3_seq = 'CACGAATAAATAATATAAGATTTTGATTTTTGCCCCCCTCATTAACTTTACTGTTAATAAGAAGAGTTGTAGAAACAGGAACTGGGACAGGTTGAACAGTTTATCCTCCTTTATCATCAATTATTGCTCATACTGGAGCCTCTGTAGATTTATCTATTTTCTCTTTTCATATAGCCGGAATTTCATCCATTCTAGGAGCTATTAACTTTATTACAACTATAATTAATATACGAGTAAAAAATATTAAATTTGACCAAATCCCTTTATTCTCCTGATCAGTAATTATTACAGCTATTTTACTTTTATTATCTCTACCAGTCTTAGCCGGAGCTATTACTATATTATTAACAGATCGAAATCTTAATACCTCTTTTTTCGACCCCAGGGGAGGAGGGGATCCTATTTTATACCAACACTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCACTCTGCTCTGACTCTCCTACC'
  t3_phred = '~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~d~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~}~~~~~~~X~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  
  test_seqdenoise3 = DNAseq(t3_seq, name = "test3", phred = t3_phred)
  test_seqdenoise3 = frame(test_seqdenoise3)
  test_seqdenoise3 = adjust(test_seqdenoise3)
  test_seqdenoise3 = outseq(test_seqdenoise3)
  
  test_seqdenoise3$name
  test_seqdenoise3$raw  
  test_seqdenoise3$data
  test_seqdenoise3$adjusted_sequence
    
})