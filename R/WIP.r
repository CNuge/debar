# # Work in progress file - for dev purposes only
# 
# # load("R/sysdata.rda")
# 

# Currently working on being able to process the data from the raw sequel data files in a way that preserves
# the phred information from input to output.
# Also trying to make the program robust to forward and backwards reads and adding logic so that it
# Tosses out bad sequences as opposed to processing them.

# Therefore 4 goals to achieve and test before actually trying full denoising of .fastq files from sujeevan.
# 1. Select the better read from forward and backwards orientations (option that can be turned off if all expected to be forward).
# Also have it so that mirror reads (i.e. around the bell of the adapter so both 5'-3' and 3'-5' in the same read) 
# are identified (by match number?) removed.
# 2. Filter input sequences based on the phred quality values. avg phred and count of x below a certain value. (options with defaults)
# 3. Carry the phred scores through along with the base pairs by making the vector of base pairs into a labelled vector where  
# the corresponding phred scores are the names.
# 4. Write a log file when the output sequence file is being generated. keep track of the number of low phred, the number of
# mixed orientiation and number of good files, possibly an option where bad reads aren't tossed but rather piped to a 'rejects file.
