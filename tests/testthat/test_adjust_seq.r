
random_bp = function(exclude_base = NULL){
  bps = c('a', 't', 'g', 'c')
  if(!is.null(exclude_base)){
    bps = bps[bps != exclude_base]
  }
  sample(bps, 1)
}


error_introduce = function(dna_string, 
                            global_mutation_rate = 0.01, 
                            global_indel_rate = 0.01){
  org_vec = tolower(strsplit(dna_string, "")[[1]])
  new_seq = c()
  point_positions = c()
  insert_positions = c()
  delete_positions = c()

  for(i in 1:length(org_vec)){
    b = org_vec[i]
    prob = runif(1)
    if(prob < global_mutation_rate){
      #point mutation
      new_b = random_bp(exclude_base=b)
      new_seq = c(new_seq,new_b)
      point_positions = c(point_positions, i)
    } else if ((global_mutation_rate < prob) && (prob<(global_indel_rate+global_mutation_rate))){
      #indel
      in_prob = runif(1)
      if(in_prob<0.5){
        #insertion
        #add the base
        new_seq = c(new_seq, b)
        #insert a base after
        new_seq = c(new_seq, random_bp())
        insert_positions = c(insert_positions, i)
      }else{
        #deletion
        #don't add anything
        delete_positions = c(delete_positions, i)
      }
    }else{
      new_seq = c(new_seq, b)
    }
  }
  return(list(
    original = dna_string,
    error = paste(new_seq, collapse = ""),
    point_positions = point_positions, 
    insert_positions = insert_positions ,
    delete_positions = delete_positions))
}


dna_string = "ctctacttgatttttggtgcatgagcaggaatagttggaatagctttaagtttactaattcgcgctgaactaggtcaacccggatctcttttaggggatgatcagatttataatgtgatcgtaaccgcccatgcctttgtaataatcttttttatggttatacctgtaataattggtggctttggcaattgacttgttcctttaataattggtgcaccagatatagcattccctcgaataaataatataagtttctggcttcttcctccttcgttcttacttctcctggcctccgcaggagtagaagctggagcaggaaccggatgaactgtatatcctcctttagcaggtaatttagcacatgctggcccctctgttgatttagccatcttttcccttcatttggccggtatctcatcaattttagcctctattaattttattacaactattattaatataaaacccccaactatttctcaatatcaaacaccattatttgtttgatctattcttatcaccactgttcttctactccttgctctccctgttcttgcagccggaattacaatattattaacagaccgcaacctcaacactacattctttgaccccgcagggggaggggacccaattctctatcaacactta"

error_introduce("ATGCA")
error_introduce(dna_string)