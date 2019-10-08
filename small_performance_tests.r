library(seqdenoise)

fastq_dat_file = system.file('extdata/ccs_subset.fastq', package = 'seqdenoise')
#original sequence data example from Hebert et al. 2018 data
fastq_dat_file_messy = system.file('extdata/previous_pacbio_example.fastq', package = 'seqdenoise')
#
#new ccs data, from 158 SMRT
fastq_real_file = system.file('extdata/sequel_smrt_subset.fastq', package = 'seqdenoise')
#new ccs data, from 158 SMRT
#in gz format
fastq_gz_real_file = system.file('extdata/sequel_smrt_subset.fastq.gz', package = 'seqdenoise')

#new ccs data, from 229 SMRT
fastq_real_file2 = system.file('extdata/sequel_smrt_subset2.fastq', package = 'seqdenoise')

#for exploriation purposes, iontorrent example
#mbrave: CBGMB-00029
fastq_ions5 = system.file('extdata/iontorrent_s5_example.fastq', package = 'seqdenoise')


denoise_file(fastq_gz_real_file, keep_flanks = TRUE, filename = "with_edges_seqdenoise_out.fastq")

denoise_file(fastq_gz_real_file, keep_flanks = FALSE, filename = "noedges_seqdenoise_out.fastq")

denoise_file(fastq_dat_file, filename = "prev_example_seqdenoise.fastq")

denoise_file(fastq_dat_file_messy, filename = "prev_example_seqdenoise2.fastq")

denoise_file(fastq_real_file2, filename = "with_edges_seqdenoise_out_example2.fastq")

denoise_file(fastq_ions5, filename = "iontorrent_example_3.fastq")

denoise_file(fastq_ions5,  keep_flanks = FALSE, filename = "iontorrent_example_3_noedges.fastq")

denoise_file(fastq_ions5,  keep_flanks = FALSE, filename = "iontorrent_example_3_noedges_with_maksed.fastq")
denoise_file(fastq_ions5, keep_flanks = FALSE, write_masked = FALSE, filename = "iontorrent_example_3_noedges_nomasked.fastq")

denoise_file(fastq_ions5, multicore=4, keep_flanks = FALSE, write_masked = FALSE, filename = "iontorrent_example_3_noedges_nomasked_multicore.fastq")


#to compress the output:
#system("gzip #THEOUTPUT NAME")

################
# explore the new outputs

new_ccs_df = read_fastq(fastq_gz_real_file)

for(i in 1:length(new_ccs_df$sequence)){

	out = seqdenoise::denoise(new_ccs_df$sequence[i], to_file= FALSE)

}


i = 2
x = coi5p_pipe(new_ccs_df$sequence[i])
x$data$ntPath
x$raw
x$framed
x

rev_ex = paste(seqinr::comp(strsplit(new_ccs_df$sequence[i], "")[[1]]), collapse="")
y = coi5p_pipe(rev_ex)
y$data$ntPath
y$raw



t2 = "gaactttgtactttatatttggagcatgagctggaatagtaggaacttctctaagaatcttaattcgagctgaattaggacatccaggagcattaattggtgatgatcaaatttataatgtaattgtcactgcccatgctttcgtaataatttttttatagtaataccaattataattggaggatttggaaattgattagttccattaatattaggagctccagatatagcatttcctcgaataaataatatgagtctctgattattaccgccttcattgactctattaatggctagaagtatagtagaaaatggggctggtactggatgaactgtttatccacctttatcttctattattgctcatggtggtgcttcagtagatttagcaattttttctttacatttagcaggagtttcttcaattttaggagcgctcctaaaattgaagaaactcctgctaaatgtaaagaaaaaattgctaaatctactgaagcaccaccatgagcaataatagaagataaaggtggataaacagttcatccagtaccagccccattttctactatacttctagccattaatagagtcaatgaaggcggtaataatcagaaactcatattatttattcgaggaaatgctatatccggagctcctaatattaatggaactaatcaatttccaaatcctccaattataattggtattactataaaaaaaattattacgaaagcatgggcagtgacaattacattataaatttgatcatcaccaattaatgctcctggatgtcctaattcagctcgaattaagattcttagagaagttcctactattccagctcatgctccaaatataaagtacaaagttccaatatctttatgattggttgactgtgacctgagtcagctacatgttcgactgcatagacagtgcgctgtctatgcagtcgaacatgtagctgactcaggtcacggtcaacaaatcataaagatattg"
x = coi5p_pipe(t2)
x$data$ntPath
x$raw
nchar(x$framed)
x
################
# test the reverse compliment hypothesis

for(i in 1:length(new_ccs_df$sequence)){

	rev_str = paste(seqinr::comp(strsplit(new_ccs_df$sequence[i], "")[[1]]), collapse="")
	seqdenoise::denoise(rev_str, name = new_ccs_df$header_data[i], 
					filename = 'rev_comp_seqdenoise_out.fastq')

}


##############
# Would the old info have failed on rev comp sequences?
# - yes 
ccs_df = read_fastq(fastq_dat_file)

for(i in 1:length(ccs_df$sequence)){

	rev_str = paste(seqinr::comp(strsplit(ccs_df$sequence[i], "")[[1]]), collapse="")

	seqdenoise::denoise(rev_str, name = ccs_df$header_data[i], 
					filename = 'rev_comp_seqdenoise_old_data_out.fastq')

}






library(coiDenoiser)

ccs_df = read_fastq(fastq_dat_file)

for(i in 1:length(ccs_df$sequence)){
	dat = coiDenoiser::build_ccs(x = ccs_df$sequence[i],
				id =ccs_df$header_data[i])

	coiDenoiser::denoise(dat, censor_length = 5, 
				to_file = TRUE, outformat = 'fastq',
				filename = 'coiDenoiser_out_previous.fastq')
}


new_ccs_df = read_fastq(fastq_gz_real_file)

for(i in 1:length(new_ccs_df$sequence)){
	dat = coiDenoiser::build_ccs(x = new_ccs_df$sequence[i],
				id =new_ccs_df$header_data[i])

	coiDenoiser::denoise(dat, censor_length = 5, to_file = TRUE, filename = 'coiDenoiser_new_data.fastq')
}
