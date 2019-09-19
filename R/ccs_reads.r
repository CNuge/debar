
#' Build a new coi5p class instance.
#'
#' @keywords internal
new_DNAseq = function(x = character(), name = character()){
  stopifnot(is.character(x))
  stopifnot(is.character(name))
  if(length(x) == 0){
    stop("Must pass a DNA sequence.")
  }
  structure(list(name = name, raw = tolower(x)) , class = "DNAseq")
}

#' Validate the new coi5p class instance.
#'
#' @keywords internal
validate_DNAseq = function(new_instance){
  # take a new instance and run validation checks on the sequence
  # make sure the sequence has only ATGCN-
  # make sure the sequence has length greater than zero
  allowed = c("-", "a", "c", "g", "n","t")
  for(c in sort(unique(strsplit(new_instance$raw, "")[[1]]))){
    if(!c %in% allowed){
      stop(paste("Unallowed character in DNA string:", c,
                 "\nValid characters are: a t g c - n"))
    }
  }
  new_instance
}


#' Build a coi5p object from a DNA sequence string.
#'
#' @param x a nucleotide string.
#' Valid characters within the nucleotide string are: a,t,g,c,-,n.
#' The nucleotide string can be input as upper case, but will be automatically converted to lower case.
#' @param name an optional character string. Identifier for the sequence.
#'
#' @return an object of class \code{"coi5p"}
#' @examples
#' dat = DNAseq(example_nt_string)
#' #named DNAseq sequence
#' dat = DNAseq(example_nt_string, name = "example_seq1")
#' #components in output DNAseq object:
#' dat$raw
#' dat$name
#' @name DNAseq
#' @export
DNAseq = function(x = character(), name = character()){
  validate_DNAseq(new_DNAseq(tolower(x), name))
}

