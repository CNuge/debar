% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjust_seq.r
\name{adj_seq}
\alias{adj_seq}
\title{Adjust the DNA sequence based on the ntPHMM path.}
\usage{
adj_seq(frame_dat, path_out, censor_length = 3, added_phred = "*")
}
\arguments{
\item{frame_dat}{The DNAseq's framing data - generated by set_frame()}

\item{path_out}{The nucleotide PHMM path for the sequence.}

\item{censor_length}{Number of base pairs in either direction of a PHMM correction
to convert to placeholder characters.}

\item{added_phred}{The phred character to use for characters inserted into the original sequence.}
}
\value{
A named list.
}
\description{
Adjust the DNA sequence based on the ntPHMM path.
}
\keyword{internal}
