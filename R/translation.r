#' Censored Translation of a codon.
#'
#' Translate a codon of DNA sequence using the censored translation table.
#' this translates codons for which the amino acids is unambigious across
#' mitochondrial genetic codes across the animal kingdom and does not
#' translate those for which the amino acid varies,
#' but rather outputs a ? in the string.
#' @param codon a three letter DNA string.
#' @details
#' Censored translation table:
#'            FFLLSSSSYY?*CCWWLLLLPPPPHHQQRRRRII?MTTTTNN?KSS??VVVVAAAADDEEGGGG
#'   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#'   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#'   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
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
#' @details
#' Censored translation table:
#'      AA  = FFLLSSSSYY?*CCWWLLLLPPPPHHQQRRRRII?MTTTTNN?KSS??VVVVAAAADDEEGGGG
#'   Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
#'   Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
#'   Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
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


#' Determine the translation table to use for a given phylogenetic group.
#'
#' Recommends which translation table to use if taxonomic data is avaliable.
#' The reccomendations are based on the translation tables reported for different
#' taxonomic classifications on the barcode of life database
#' (http://www.boldsystems.org/index.php).
#'
#' @param x a taxonomic designation (allowed ranks: family, order, class, phylum).
#' @return an integer indicating the correct translation table.
#' @examples
#' which_trans_table("Chordata") #phylum
#' which_trans_table("Actinopterygii") #class
#' which_trans_table("Akentrogonida")  #order
#' which_trans_table("Hydrobiidae") #family
#' @details
#' If which_trans_table is unable to identify a translation table to utilize,
#' more information on translation tables can be found here: https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
#' @export
which_trans_table = function(x) {
    use_tab = trans_df$trans_table[trans_df$phylogeny == tolower(x)]
    if(length(use_tab) == 0){
      return(0)
    }
    return(use_tab)
}
