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





#TODO - need to stop masking for triple deletes, these are valid missing codons.
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

