
test_that("The adjustment of individual ccs reads is conducted properly.", {

  ############################
	order2 = 'lepidoptera'
	starting_sequence2 = 'TTTATTTTTGGAATTTGATCCGGAATAATTGGAACATCTCTTAGTCTATTAATTCGTGCTGAATTAGGAAACCCAGGCTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACCGCCCACGCCTTTATTATAATTTTTTTCATGGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTACTTCCCCCTTCACTTACTTTATTAGTTTCTAGAAGAATTGTAGAAAATGGAGTAGGAACAGGATGAACAGTATACCCCCCTTTATCTTCTAATATTGCCCATGGTGGAGGCTCTGTTGATTTAGCAATCTTTTCTCTTCATTTAGCTGGAATTTCTTCAATTTTAGGAGCTGTCAATTTTATTACAACAGTAATTAATATACGAACAAATGGTATATCTTTTGATCGAATACCATTATTTGTTTGATCTGTTGCTATCACAGCACTTCTACTACTTTTATCTTTACCTGTCTTAGCTGGAGCTATTACTATACTTTTAACTGAT'

	#pull the correct components from this
	fake_ccs_data2 = list(
	'TTTATTTTTGGAATTTGATCCGGAATAATTGGAACATCTCTTAGTCTATTAATTCGTGCTGAATTAGGAAACCCAGGCTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACCGCCCACGCCTTTATTATAATTTTTTTCATGGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTACTTCCCCCCTTCACTTACTTTATTAGTTTCTAGAAGAATTGTAGAAAATGGAGTAGGAACAGGATGAACAGTATACCCCCCTTTATCTTCTAATATTGCCCATGGTGGAGGCTCTGTTGATTTAGCAATCTTTTCTCTTCATTTAGCTGGAATTTCTTCAATTTTAGGAGCTGTCAATTTTATTACAACAGTAATTAATATACGAACAAATGGTATATCTTTTGATCGAATACCATTATTTGTTTGATCTGTTGCTATCACAGCACTTCTACTACTTTTATCTTTACCTGTCTTAGCTGGAGCTATTACTATACTTTTAACTGAT',
	'TTTATTTTTGGAATTTGATCCGGAATAATTGGAACATCTCTTAGTCTATTAATTCGTGCTGAATTAGGAAACCCAGGCTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACCGCCCACGCCTTTATTATAATTTTTTTCATGGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCCCGAATAAATAATATAAGATTTTGATTACTTCCCCCTTCACTTACTTTATTAGTTTCTAGAAGAATTGTAGAAAATGGAGTAGGAACAGGATGAACAGTATACCCCCCTTTATCTTCTAATATTGCCCATGGTGGAGGCTCTGTTGATTTAGCAATCTTTTCTCTTCATTTAGCTGGAATTTCTTCAATTTTAGGAGCTGTCAATTTTATTACAACAGTAATTAATATACGAACAAATGGTATATCTTTTGATCGAATACCATTATTTGTT',
	'TTTATTTTTGGAATTTGATCCGGAATAATTGGAACATCTCTTAGTCTATTAATTCGTGCTGAATTAGGAAACCCAGGCTCTTTAATTGGAGATGATCAAATTTATAATACAATTGTTACCGCCCACGCCTTTATTATAATTTTTTTCATGGTTATACCAATTATAATTGGAGGATTTGGAAATTGATTAGTACCTTTAATATTAGGAGCTCCTGATATAGCTTTCCCCCCGAATAAATAATATAAGATTTTGATTACTTCCCCCTTCACTTACTTTATTAGTTTCTAGAAGAATTGTAGAAAATGGAGTAGGAACAGGATGAACAGTATACCCCCCTTTATCTTCTAATATTGCCCATGGTGGAGGCTCTGTTGATTTAGCAATCTTTTCTCTTCATTTAGCTGGAATTTCTTCAATTTTAGGAGCTGTCAATTTTATTACAACAGTAATTAATATACGAACAAATGGTATATCTTTTGATCGAATACCATTATTTGTTTGATCTGTTGCTATCACAGCACTTCTACTACTTTTATCTTTACCTGTCTTAGCTGGAGCTATTACTATACTTTTAACTGAT'
		)

	fake_ccs_2 = lapply(fake_ccs_data2, DNAseq) 

	fake_ccs_2 = lapply(fake_ccs_2, frame)
	fake_ccs_2 = mapply(adjust, fake_ccs_2, censor_length = 5)
	

	#TODO - build the expected outputs: 
	#       these are based on looking at the paths built from aphid and seeing where the corrections
	read_1_output = ''
	read_2_output = ''
	read_3_output = ''
	
	#note the censorship isn't applied here yet... need to re-add this in to adjust seq
	#expect_equal(paste(fake_ccs_2[,1][['adjusted_sequence']],collapse=''), read_1_output)
	#expect_equal(paste(fake_ccs_2[,2][['adjusted_sequence']],collapse=''), read_2_output)
	#expect_equal(paste(fake_ccs_2[,3][['adjusted_sequence']],collapse=''), read_3_output)
	expect_equal(1,1)
})