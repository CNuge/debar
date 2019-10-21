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
  
  
  ########################
  # This one has low QV and should therefore be rejected
  t4_seq = "AGTTAACAAATCATAAAGATATTGGTACCTTATACTTTATTTTTGGAGCTTGAGCTGGTATAATTGGCACATCCCTTAGTATTATTATTCGTGCAGAATTAGGACATCCAGGTACATTTATTGGTAATGACNNNNNNNNNNNTGTTATTGTTACTGCTCATGCCTTTATTATAATTTTTTTCATAGTTATACCAATTATAATCGGTGGATTTGGAAATTGACTTGTACCTTTAATATTAGGTGCCCCCGATATAGCCTTTCCCCNNNNNNNNNNTATAAGATTTTGAATATTACCCCCTTCTATTACTTTATTACTAACAAGAAGACTTGTAGAGAATGGAGCAGGTACCGGATGAACAGTTTACCCACCATTAGCCGCTAATATTTCTCATGCTGGATCCTCAGTTGACCTAGCAATTTTCTCTTTACACTTAGCTGGGATTTCTTCCATTCTAGGAGCTGTAAATTTTATTACTACAATTATCAATATACGATCAAATGGAATTTCTTTCGACCGTATGCCTTTATTTGTATGATCTGTCTTTATTACAGCTATTTTATTACTTTTATCATTGCCTGTTTTAGCCGGAGCAATTACAATATTATTAACTGATCGTAATTTAAATACATCATTCTTTGACCCTGCAGGAGGGGGAGATCCAATTCTTTACCAACATTTATTTTGATTTTTGGTCATCCATGAAGTTTA"
  t4_phred = "f<$)m)a.h4fvhfn{;mw[]z^z^zeyqz$f[anxi@wjpu$1fycjn~pmtsjzjqlxd_@tndylkbmxw{w|ufkczea$vr[g&lc9>aa_i+chb~=$qgnzbc|[gfoc]m4d|]~$kx@ar$pd=q+*l0c_r0$j~ke}]johyykhsl^[9d|aem*ymms_ci[x|{wr)]ep^g}anfgqyglr_epojypjzjmyizonitr~?z$qap}[^ri|iap:jnkf~n_~]h_przt3_]_zbe^sp$t3syx1ockj-l-jg>awmxd^i|||4@>gkebc]qfkx7$aeuq_9f[cx+>w[ng$qcyerornojy^:ebaqgvigh}q$nsaie]wkmlb7vmo=xfh{{4ezz@g$2>_>gpduqv$d8aslz=_kreydpogohj}shpfbx;gh_=i=l^pzu{|>[jppudcrmozuyypviked_h>oe5$0nw^qbd_&$axh_jdt0[us?btbdqbg[c9yzxeqrbtpk$^nysxcjau]c{zyt9jijefrjthkzsx[vnce;cqi<n-w>$]*x<lnl]q99jopiqhjas[qq2bvdpitly/gp]d|tni:winz{:ci}rw[hm[ub_ktmfw_evvdxwlouipinqjn_lisfsb)jcvrf]dhl_^$:w[]gfjb$wx?als{{1cj$kuxlhszp_qdr[:gex[la7$d_a$nbkxws0$.$5kvb;-plhaclh>^"
  
  test_seqdenoise4 = DNAseq(t4_seq, name = "test4", phred = t4_phred)
  
  test_seqdenoise4 = frame(test_seqdenoise4)
  test_seqdenoise4 = adjust(test_seqdenoise4)
  test_seqdenoise4 = outseq(test_seqdenoise4)
  
  test_seqdenoise4$name
  test_seqdenoise4$raw  
  test_seqdenoise4$data
  test_seqdenoise4$adjusted_sequence
  
  
  #####################
  # This one appears to be dropping the phred scores in the fq output, see if you can pin down the bug
  @m54146_190221_235805/15859816/ccs
  t5_seq = "GAGTCAGAGATCTATGGATCACTTGTGCAAGCATCACATCGTAGTAAACTTCTGGGTGTCCAAAAAAATCAAAATAAATGTTGATATAAAATTGGATCACCTCCTCCTGCTGGGTCAAAGAAAGATGTATTTAAATTTCGGTCGGTTAATAATATTGTAATAGCTCCAGCTAATACAGGTAGAGATAAAAGTAATAAAATAGCAGTAATTACAACAGATCATACAAATAAAGGTATTCGATCAAATGTAATACCTGTTGATCGTATATTAATTACAGTAGTAATGAAATTTACGGCTCCTAAAATAGAAGAAATACCAGCTAAATGTAATGAAAAAAATAGCTAAATCAACTGAGGCTCCTCCATGAGCAATCCCTGCTGATAGAGGAGGGTAAACTGTTCAACCAGTTCCAGCTCCATTTTCTACTATACTACTAGCTAATAATAATGTAAGAGAAGGGGGTAATATTCAAAAACTTATATTACCCCCTTCTCTTACATTATTATTAGCTAGTAGTATAGTAGAAAATGGAGCTGGAACTGGTTGAACAGTTTACCCTCCTCTATCAGCAGGGATTGCTCATGGAGGAGCCTCAGTTGATTTAGCTATTTTTTCATTACATTTAGCTGGTACTTCTTCTATTTTAGGAGCCGTAAATTTCATTACTACTGTAATTAATATACGATCAACAGGTATTACATTTGATCGAATACCTTTATTTGTATGATCTGTTGTAATTACTGCTATTTTATTACTTTTATCTCTACCTGTATTAGCTGGAGCTATTACAATATTATTAACCGACCGAAATTTAAATACATCTTTCTTTGACCCAGCAGGAGGAGGTGATCCAATTTTATATCAACATTTATTTTGATTTTTTGGACATCCAGAAGTTTACTACGATGTGATGCTTGCACAAGTGATCCATAGATCTCTGACTC"
  t5_phred = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~V~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~0~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~t~~~~~~~~~~~~a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~V~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

  test_seqdenoise5 = DNAseq(t5_seq, name = "test5", phred = t5_phred)
  
  test_seqdenoise5 = frame(test_seqdenoise5)
  test_seqdenoise5 = adjust(test_seqdenoise5)
  test_seqdenoise5 = outseq(test_seqdenoise5)
  
  test_seqdenoise5$name
  test_seqdenoise5$raw  
  test_seqdenoise5$data
  test_seqdenoise5$adjusted_sequence
  
  
})