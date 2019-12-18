#' Censored Translation of a codon.
#'
#' Translate a codon of DNA sequence using the censored translation table.
#' this translates codons for which the amino acids is unambigious across
#' mitochondrial genetic codes across the animal kingdom and does not
#' translate those for which the amino acid varies,
#' but rather outputs a ? in the string.
#' @param codon a three letter DNA string.
#' @return a string.
#' @keywords internal
translate_codon = function(codon){
	trans_table = list(
		'TTT' = 'F', 'TTC' = 'F', 'TTA' = 'L', 'TTG' = 'L',
		'TCT' = 'S', 'TCC' = 'S', 'TCA' = 'S', 'TCG' = 'S',
		'TAT' = 'Y', 'TAC' = 'Y', 'TAA' = '?', 'TAG' = '*',
		'TGT' = 'C', 'TGC' = 'C', 'TGA' = 'W', 'TGG' = 'W',
		'CTT' = 'L', 'CTC' = 'L', 'CTA' = 'L', 'CTG' = 'L',
		'CCT' = 'P', 'CCC' = 'P', 'CCA' = 'P', 'CCG' = 'P',
		'CAT' = 'H', 'CAC' = 'H', 'CAA' = 'Q', 'CAG' = 'Q',
		'CGT' = 'R', 'CGC' = 'R', 'CGA' = 'R', 'CGG' = 'R',
		'ATT' = 'I', 'ATC' = 'I', 'ATA' = '?', 'ATG' = 'M',
		'ACT' = 'T', 'ACC' = 'T', 'ACA' = 'T', 'ACG' = 'T',
		'AAT' = 'N', 'AAC' = 'N', 'AAA' = '?', 'AAG' = 'K',
		'AGT' = 'S', 'AGC' = 'S', 'AGA' = '?', 'AGG' = '?',
		'GTT' = 'V', 'GTC' = 'V', 'GTA' = 'V', 'GTG' = 'V',
		'GCT' = 'A', 'GCC' = 'A', 'GCA' = 'A', 'GCG' = 'A',
		'GAT' = 'D', 'GAC' = 'D', 'GAA' = 'E', 'GAG' = 'E',
		'GGT' = 'G', 'GGC' = 'G', 'GGA' = 'G', 'GGG' = 'G'
		)

	if(nchar(codon) != 3){
		return('-')
	}

	if(grepl('-', codon,  fixed = TRUE)){
		return('-')
	}

	if(grepl('N', (codon),  fixed = TRUE)){
		return('-')
	}

	if(grepl('n', (codon),  fixed = TRUE)){
	  return('-')
	}
	return(trans_table[[toupper(codon)]])
}

#' Censored Translation of a DNA string.
#'
#' Translate a DNA sequence using the censored translation table,
#' this translates codons for which the amino acids is unambigious across
#' mitochondrial genetic codes across the animal kingdom and does not
#' translate those for which the amino acid varies,
#' but rather outputs a ? in the string.
#' @param dna_str The DNA string to be translated.
#' @param reading_frame reading frame = 1 means the first bp in the string is the start of the
#' first codon, can pass 1, 2 or 3. For 2 and 3 the first 1 and 2 bp will be
#' dropped from translation respectively.
#' @seealso \code{\link{aa_check}}
#' @export
censored_translation = function(dna_str, reading_frame = 1){
	num_bp = nchar(dna_str)

	codons = seq(reading_frame, num_bp, by=3)

	codon_vec = sapply(codons, function(x) {
						substr(dna_str, x, x+2)
						})

	aa_str = paste(lapply(codon_vec, translate_codon), collapse= "")

	return(aa_str)
}

