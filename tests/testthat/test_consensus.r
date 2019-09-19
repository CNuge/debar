test_that("A consensus sequence is obtained properly", {


	expect_equal(col_mode(c(NA, NA, NA)),'n')
	expect_equal(col_mode(c(NA, NA, NA, 't')) ,'t')
	expect_equal(col_mode(c(NA, NA, 't', 't')),'t')
	expect_equal(col_mode(c(NA, 't', 't', 't')),'t')
	expect_equal(col_mode(c(NA, 't', 't', 'g')) ,'t')
	expect_equal(col_mode(c(NA, 't', 't', 'c', 'g')), 't')
	expect_equal(col_mode(c(NA, 't', 't', 'g','g')),'n')
	expect_equal(col_mode(c(NA, 't', 't' ,'g' ,'g', 'c', 'c')) ,'n')


	expected_consensus_sequence = "NNNNNNNNNTTTATTTTTGGAATTTGATCCGGAATAATTGGAACATCTCTTAGTCTATTAATTCGTGCTGAATTAGGAAACCCAGGCTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACCGCCCACGCCTTTATTATAATTTTTTTCATGGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTACTTCCCCCTTCACTTACTTTATTAGTTTCTAGAAGAATTGTAGAAAATGGAGTAGGAACAGGATGAACAGTATACCCCCCTTTATCTTCTAATATTGCCCATGGTGGAGGCTCTGTTGATTTAGCAATCTTTTCTCTTCATTTAGCTGGAATTTCTTCAATTTTAGGAGCTGTCAATTTTATTACAACAGTAATTAATATACGAACAAATGGTATATCTTTTGATCGAATACCATTATTTGTTTGATCTGTTGCTATCACAGCACTTCTACTACTTTTATCTTTACCTGTCTTAGCTGGAGCTATTACTATACTTTTAA"

	list_of_adjusted_seqs = list(
		unlist(strsplit("----------------------------------------gaacatctcttagtctattaattcgtg------taggaaacccaggctctttaattggagatgatcaaatttataatacaattgttaccgcccacgcctttattataatttttttcatggttataccaattataattggaggatttggaaattgattagtacctttaatattaggagctcctgatatagctttcccccgaataaataatataagattttgattacttcccccttcacttactttattagtttctagaagaattgtagaaaatggagtaggaacaggatgaacagtatacccccctttatcttctaatattgcccatggtggaggctctgttgatttagcaatcttttctcttcatttagctggaatttcttcaattttaggagctgtcaattttattacaacagtaattaatatacgaacaaatggtatatcttttgatcgaataccattatttgtttgatctgttgctatcacagcacttctactacttttatctttacctgtcttagctggagctattactatacttttaa", split = "")),
		unlist(strsplit("---------tttatttttggaatttgatccggaataattggaacatctcttagtctattaattcgtgctgaattaggaaacccaggctctttaattggagatgatcaaatttataatacaattgttaccgcccacgcctttattataatttttttcatggttataccaattataattggaggatttggaaattgattagtacctttaatattaggagctcctgat------ttcccccgaataaataatataagattttgattacttcccccttcacttactttattagtttctagaagaattgtagaaaatggagtaggaacaggatgaacagtatacccccctttatcttctaatattgcccatggtggaggctctgttgatttagcaatcttttctcttcatttagctggaatttcttcaattttaggagctgtcaattttattacaacagtaattaatatacgaacaaatggtatatcttttgatcgaataccattatttgtttgatctgttgctatcacagcacttctactacttttatctttacctgtcttagctggagctattactatacttttaa", split = "")),
		unlist(strsplit("---------tttatttttggaatttgatccggaataattggaacatctcttagtctattaattcgtgctgaattaggaaacccaggctctttaattggagatgatcaaatttataatacaattgttaccgcccacgcctttattataatttttttcatggttataccaattataattggaggatttggaaattgattagtacctttaatattaggagctcctgatatagctttcccccgaataaataatataagattttgattacttcccccttcacttactttattagtttctagaagaattgtagaaaatggagtaggaacaggatgaacagtatacccccctttatcttctaatattgcccatggtggaggctctgttgatttagcaatcttttctcttcatttagctggaatttcttcaattttaggagctgtcaattttattacaacagtaattaatatacgaacaaatggtatatcttttga", split = "")),
		unlist(strsplit("---------tttatttttggaatttgatccggaataattggaacatctcttagtctattaattcgtgctgaattaggaaacccaggctctttaattggagatgatcaaatttataatacaattgttaccgcccacgcctttattataatttttttcatggttataccaattataattggaggatttggaaattgattagtacctttaatattaggagctcctgatatagctttcccccgaataaataatataagattttgattacttcccccttcacttactttattagtttctagaagaattgtagaaaatggagtaggaacaggatgaacagtatacccccctttatcttctaatattgcccatggtggaggctctgt------agcaattttttctcttcatttagctggaatttcttcaattttaggagctgtcaattttattacaacagtaattaatatacgaacaaatggtatatcttttgatcgaataccattatttgtttgatctgttgctatcacagcacttctact", split = ""))
		)

	 #consensus_sequence(list_of_adjusted_seqs) == expected_consensus_sequence
	expect_equal(consensus_sequence(list_of_adjusted_seqs),  expected_consensus_sequence)

})
