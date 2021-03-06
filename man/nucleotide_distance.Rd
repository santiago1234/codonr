% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_codons.R
\name{nucleotide_distance}
\alias{nucleotide_distance}
\title{Nucleotide distance}
\usage{
nucleotide_distance(seq1, seq2)
}
\arguments{
\item{seq1}{string dna sequence}

\item{seq2}{dna sequence to compare same length as seq1}
}
\value{
integer number of different characters
}
\description{
How many nucleotides (characters) are different between the two strings?
}
\examples{
nucleotide_distance("AA", "AC")
}
