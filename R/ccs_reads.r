
#' Build ccs_reads class instance.
#'
#' @keywords internal
new_ccs_reads = function(x = c(), id = character(), order= "UNK"){
  stopifnot(is.vector(x))
  stopifnot(is.character(id))
  stopifnot(is.character(order))
  if(length(x) == 0){
    stop("Missing argument. Must pass a DNA sequence, or vector of sequences.")
  }
  
  structure(list(id = id, 
                 sequences = x,
                 order = tolower(order)) , class = "ccs_reads")
}

#' Validate the new ccs_reads class instance.
#'
#' @keywords internal
validate_ccs_reads = function(new_instance){
  #add sanity checks to make sure there aren't errors in the ccs_read strucutre
  #that was just built
  new_instance
}

#' Build a ccs_reads object from a DNA sequence string.
#'
#' @param x a vector of nucleotide strings.
#' @param id and optional identification string.
#' @param order the taxonomic order corresponding to the sample. 
#' If taxonomy is not known, use the default 'unknown' option.
#' @return an object of class code{"coi5p"}
#' @examples
#' #known taxonomy:
#' ex_data = build_ccs(ex_ccs_read_list, order = 'Diptera', id = 'custom_id')
#' #unknown taxonomy:
#' ex_data = build_ccs(ex_ccs_read_list,  id = 'custom_id')
#' @name build_ccs
#' @export
build_ccs = function(x, id = character(), order = "unknown"){
  validate_ccs_reads(new_ccs_reads(x, id = id, order = order))
}
