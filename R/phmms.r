#load('R/sysdata.rda')
#use_data(nt_PHMM, example_nt_string, example_nt_string_errors,ex_nt_list, overwrite = TRUE, internal = TRUE)
#then load for dev with:
# load('R/sysdata.rda')


###############################################################################
#'An example of a list of four coi5p sequences, each containing indel errors.
#'
#'
"ex_nt_list"
###############################################################################

###############################################################################
#' Nucleotide profile hidden markov model for debar.
#'
#' This model is stored in the package and was trained on a representitive
#' sample of the barcode of life database (http://www.boldsystems.org/index.php).
#'
#' @keywords internal
"nt_PHMM"
###############################################################################

###############################################################################
#' Example coi5p DNA sequence string.
#'
#' This string of barcode data is used in the package documentation's examples
#' and within the vignette demonstrating how to use the package.
"example_nt_string"
###############################################################################


###############################################################################
#' Example coi5p DNA sequence string with insertion and deletion errors.
#'
#' This string of barcode data is used in the package documentation's examples
#' and within the vignette demonstrating how to use the package.
"example_nt_string_errors"
###############################################################################
