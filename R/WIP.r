# Work in progress file - for dev purposes only


# Notes on changes to make based on sujeevan's feedback
# Possibly bifurcate into a read by read code base
# 1. rejig things so that the objects are single reads as opposed to groups of reads
# 2. don't use the order by order models, instead use just a single overall model for corrections.
#    this is to prevent overfit
# 3. don't need to take the consensus of the grouped reads, this can be done later in the MBRAVE pipeline
#    may need to do this later...
# 4. retain the edges of the sequences. make the corrections in the middle but then paste the edges back on in order to 
#    create an output that retains the primer information and other adjacent info from the original sequence.
#    Hopefully this can buff out from of the errors in the middle of the sequences and that the MBRAVE outputs
#    will be more correct.



#this is branched off from coiDenoiser
# a minimal use case of only denoising one sequence at a time - no need to take the consensus 
# or to act on multiple sequences individually.

#changing the object so that it is just a single DNA sequence
#using the small PHMM (10K training) that is generic for all the different data... can sub in the higher trained 
#order specific later.



#DONE - need to stop masking for triple deletes, these are valid missing codons.
#TODO -  adjust seq seems to be dropping the final bp... is this because its being lost in the frame step?
change in adjust seq from:
  #this prevents additions beyond the end of the original sequence.
  if(org_seq_pos >= length(org_seq_vec)){
    break
  }

to:
  #this prevents additions beyond the end of the original sequence.
  if(org_seq_pos > length(org_seq_vec)){
    break
  }
#^fix this in coiDenoiser as well.

#TODO - don't need to add the information at the start of the adjusted sequence.
#TODO - use frame to just store the start and end positions in three bits:
#  $front
#  $framed
#  $end
# don't want to toss the excess information, (may need to change the double frame logic)
# instead we want this put to the side so it can be patched back on the sequence after the middle
# is denoised
# TODO - AA check with the PHMM, project the masks back down a level corresponding to the first 0 or 2 seen.
# TODO - rejig so that instead of making additional strings, just saving index positions on the initial raw sequence.
# maybe will be faster if we don't touch the characters repeatedly but instead just manipulate a vector via index positions.
# TODO - less important than making it work but will have to rewrite the documentation to account for the new data
# structures being used here
# TODO - sujeevan said that some of the reads may be the reverse compliments? figure out how these are handled... can it be fixed?



#TODO - move the new adj_seq  + triple_in functions into coiDenoiser... add these uit tests as well.
path_test1 = c(1,1,2,1,1,1,1,1,0,0,0,1,0,1,1,1,1) #return c(9,10,11)
path_test2 = c(1,1,2,1,1,1,1,1,0,0,1,0,1,1,1,1) #return null
path_test3 = c(1,1,2,1,1,1,1,1,0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,1,1,1,1,1,0,0,0,0) #return c(9,10,11,21,22,23)

triple_ins(path_test1, 1, length(path_test1))
triple_ins(path_test2, 1,length(path_test2))
triple_ins(path_test3, 1,(length(path_test3)-5))



x= DNAseq(example_nt_string)
x
x = frame(x)
x$data
x = adjust(x, censor_length = 0)
x$adjusted_sequence





err = "ctctacttgatttttggtgcatgagcaggaatagttggaatagctttaagtttactaattcgcgctgaactaggtcaacccggatctcttttagggggatgatcagatttataatgtgatcgtaaccgcccatgcctttgtaataatcttttttatggttatacctgtaataattggtggctttggcaattgacttgttcctttaataattggtgcaccagatatagcattccctcgaataaataatataagtttctggcttcttcctccttcgttcttacttctcctggcctccgcaggagtagaagctggagcaggaaccggatgaactgtatatcctcctttagcaggtaatttagcacatgctggcccctctgttgatttagccatcttttcccttcatttggccggtatctcatcaattttagcctctattaattttattacaactattattaatataaaacccccaactatttctcaatatcaaacaccattatttgtttgatctattcttatcaccactgttcttctactccttgctctccctgttcttgcagccggaattacaatattattaacagaccgcaacctcaacactacattctttgaccccgcagggggaggggacccaattctctatcaacactta"

y= DNAseq(err)
y
y = frame(y)
y$data
y = adjust(y, censor_length = 3)
y$adjusted_sequence #looks good but the last character is dropped!




