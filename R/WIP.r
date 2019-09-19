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




