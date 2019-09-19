#Use the following to create the ccs objects for unittests and make sure that the outputs
#of the pipleline steps (with the exeption of the print part) lead to the proper consensus
#sequences

#tests on real Pacbio CCS outputs from the file:
#ccs_tibble = read_tsv('/home/cnuge/bin/DAPR/data/pacbio-ccs-test-run/PACBB9380_9999_sequence_pairs.tsv')
# library(tidyverse)
ccs_tibble = read_tsv('PACBB9380_9999_sequence_pairs.tsv')

ccs_df = data.frame(ccs_ids = ccs_tibble$pb_process_id,
	sequence = ccs_tibble$pb_sequence,
	order = ccs_tibble$order,
	stringsAsFactors = FALSE  
	)

#head(ccs_df)

#split the dataframe into a list of dfs
list_of_dfs = split(ccs_df, ccs_df$ccs_ids)


ccs_list = lapply(1:length(list_of_dfs), function(i){
	build_ccs(x =list_of_dfs[[i]][["sequence"]],
						id = names(list_of_dfs[i]), 
						order = list_of_dfs[[i]][1, "order"])	
	})

# TODO -subset out the following examples and remove the call to read_tsv above

weird1_raw = list_of_dfs[['PACBB002-16']]

weird1_ccs	= build_ccs(x =weird1_raw[["sequence"]],
						id = "PACBB002-16", 
						order = weird1_raw[1, "order"])	


x = weird1_ccs
ambig_char = "N"
censor_length = 3
aa_check = FALSE
outformat = "fasta" 
filename = "denoised/assessing_singles.fasta"


weird1_ccs = denoise(weird1_ccs, 
				ambig_char = "N",
		 		censor_length = 3,
		 		aa_check = FALSE, 
		 		outformat = "fasta", 
		 		filename = "denoised/assessing_singles.fasta")

weird1_ccs

#use these to pull out the problem parts.
weird1_ccs$sequence[1]
weird1_ccs$sequence[2]

paste(weird1_ccs$data$ntPHMMout[[1]][['path']],collapse='')
paste(weird1_ccs$data$ntPHMMout[[2]][['path']],collapse='')


paste(weird1_ccs$adjusted_sequences[[1]], collapse="") 
paste(weird1_ccs$adjusted_sequences[[2]], collapse="") 


#there are only two reads
# PACBB007-16

weird2_raw = list_of_dfs[['PACBB007-16']]

weird2_ccs	= build_ccs(x =weird2_raw[["sequence"]],
						id = "PACBB002-16", 
						order = weird2_raw[1, "order"])	

weird2_ccs = denoise(weird2_ccs, 
				ambig_char = "N",
		 		censor_length = 3,
		 		aa_check = FALSE, 
		 		outformat = "fasta", 
		 		filename = "denoised/assessing_singles.fasta")

weird2_ccs

#use these to pull out the problem parts.
weird2_ccs$sequence[1]
weird2_ccs$sequence[2]
weird2_ccs$sequence[3]
weird2_ccs$sequence[4]

paste(weird2_ccs$data$ntPHMMout[[1]][['path']],collapse='')
paste(weird2_ccs$data$ntPHMMout[[2]][['path']],collapse='')
paste(weird2_ccs$data$ntPHMMout[[3]][['path']],collapse='')
paste(weird2_ccs$data$ntPHMMout[[4]][['path']],collapse='')


paste(weird2_ccs$adjusted_sequences[[1]], collapse="") 
paste(weird2_ccs$adjusted_sequences[[2]], collapse="") 
paste(weird2_ccs$adjusted_sequences[[3]], collapse="") 
paste(weird2_ccs$adjusted_sequences[[4]], collapse="") 




#'PACBB8695-16'
# this one has three codons missing relative to the profile?
# should this be accepted? seems unlikely..
weird3_raw = list_of_dfs[['PACBB8695-16']]

weird3_ccs	= build_ccs(x =weird3_raw[["sequence"]],
						id = "PACBB002-16", 
						order = weird3_raw[1, "order"])	

weird3_ccs = denoise(weird3_ccs, 
				ambig_char = "N",
		 		censor_length = 3,
		 		aa_check = FALSE, 
		 		outformat = "fasta", 
		 		filename = "denoised/assessing_singles.fasta")

weird3_ccs

#use these to pull out the problem parts.
weird3_ccs$sequence[1]
weird3_ccs$sequence[2]
weird3_ccs$sequence[3]

paste(weird3_ccs$data$ntPHMMout[[1]][['path']],collapse='')
paste(weird3_ccs$data$ntPHMMout[[2]][['path']],collapse='')
paste(weird3_ccs$data$ntPHMMout[[3]][['path']],collapse='')

paste(weird3_ccs$adjusted_sequences[[1]], collapse="") 
paste(weird3_ccs$adjusted_sequences[[2]], collapse="") 
paste(weird3_ccs$adjusted_sequences[[3]], collapse="") 


#this one for comparison with in the visualizing script
weird4_raw = list_of_dfs[['PACBB051-16']]

weird4_ccs	= build_ccs(x =weird4_raw[["sequence"]],
						id = "PACBB002-16", 
						order = weird4_raw[1, "order"])	

weird4_ccs = denoise(weird4_ccs, 
				ambig_char = "N",
		 		censor_length = 3,
		 		aa_check = FALSE, 
		 		outformat = "fasta", 
		 		filename = "denoised/assessing_singles.fasta")

weird4_ccs

#use these to pull out the problem parts.
weird4_ccs$sequence[1]
weird4_ccs$sequence[2]
weird4_ccs$sequence[3]

paste(weird4_ccs$data$ntPHMMout[[1]][['path']],collapse='')
paste(weird4_ccs$data$ntPHMMout[[2]][['path']],collapse='')
paste(weird4_ccs$data$ntPHMMout[[3]][['path']],collapse='')
paste(weird4_ccs$data$ntPHMMout[[4]][['path']],collapse='')

paste(weird4_ccs$adjusted_sequences[[1]], collapse="") 
paste(weird4_ccs$adjusted_sequences[[2]], collapse="") 
paste(weird4_ccs$adjusted_sequences[[3]], collapse="") 
paste(weird4_ccs$adjusted_sequences[[4]], collapse="") 

