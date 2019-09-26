#load('R/sysdata.rda')

#TODO - pick the  final public data model to save
# nt_filename = '/home/cnuge/bin/DAPR/public_only_PHMMs/nt_train_all_bold_coi_v2_1perc_subsample-PublicOnly.PHMM'
# new_nt = readPHMM(nt_filename)
# nt_PHMM = new_nt

# aa_filename = '/home/cnuge/bin/DAPR/public_only_PHMMs/aa_train_all_bold_coi_v2_1perc_subsample-PublicOnly.PHMM'
# new_aa = readPHMM(aa_filename)
# aa_PHMM = new_aa

#use_data(nt_PHMM, aa_PHMM, trans_df, example_nt_string, example_nt_string_errors, overwrite = TRUE, internal = TRUE)

#then load for dev with:
# load('R/sysdata.rda')



###############################################################################
#' Nucleotide profile hidden markov model for seqdenoise.
#'
#' This model is stored in the coi5p package and was trained on a representitive
#' sample of the barcode of life database (http://www.boldsystems.org/index.php).
#'
#' @keywords internal
"nt_PHMM"
###############################################################################

###############################################################################
#' Amino acid profile hidden markov model for seqdenoise.
#'
#' This model is stored in the coi5p package and was trained on a representitive
#' sample of the barcode of life database (http://www.boldsystems.org/index.php).
#'
#' @keywords internal
"aa_PHMM"
###############################################################################

###############################################################################
#' Data frame containing the translation table recommendations.
#'
"trans_df"
###############################################################################

###############################################################################
#' Example coi5p DNA sequence string
#'
#' This string of barcode data is used in the package documentation's examples
#' and within the vignette demonstrating how to use the package.
#' The variable is:
#' example_nt_string = 'ctctacttgatttttggtgcatgagcaggaatagttggaatagctttaagtttactaattcgcgctgaactaggtcaacccggatctcttttaggggatgatcagatttataatgtgatcgtaaccgcccatgcctttgtaataatcttttttatggttatacctgtaataattggtggctttggcaattgacttgttcctttaataattggtgcaccagatatagcattccctcgaataaataatataagtttctggcttcttcctccttcgttcttacttctcctggcctccgcaggagtagaagctggagcaggaaccggatgaactgtatatcctcctttagcaggtaatttagcacatgctggcccctctgttgatttagccatcttttcccttcatttggccggtatctcatcaattttagcctctattaattttattacaactattattaatataaaacccccaactatttctcaatatcaaacaccattatttgtttgatctattcttatcaccactgttcttctactccttgctctccctgttcttgcagccggaattacaatattattaacagaccgcaacctcaacactacattctttgaccccgcagggggaggggacccaattctctatcaacactta'
"example_nt_string"
###############################################################################


###############################################################################
#' Example coi5p DNA sequence string with insertion and deletion errors.
#'
#' This string of barcode data is used in the package documentation's examples
#' and within the vignette demonstrating how to use the package.
#' The variable is:
#' example_nt_string_errors = 'ctctacttgatttttggtgcatgagcaggaatagttggaatagctttaagtttactaattcgcgctgaactaggtcaaccggatctcttttaggggatgatcagatttataatgtgatcgtaaccgcccatgcctttgtaataatcttttttatggttatacctgtaataattggtggctttggcaattgacttgttcctttaataattggtgcaccagatatagcattccctcgaataaataatataagtttctggcttcttcctccttcgttcttacttctcctggcctccgcaggagtagaagctggagcaggaaccggatgaactgtatatcctccttttagcaggtaatttagcacatgctggcccctctgttgatttagccatcttttcccttcatttggccggtatctcatcaattttagcctctattaattttattacaactattattaatataaaacccccaactatttctcaatatcaaacaccattatttgtttgatctattcttatcaccactgttcttctactccttgctctccctgttcttgcagccggaattacaatattattaacagaccgcaacctcaacactacattctttgaccccgcagggggaggggacccaattctctatcaacactta'
"example_nt_string_errors"
###############################################################################
