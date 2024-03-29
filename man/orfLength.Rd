% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orfs.R
\name{orfLength}
\alias{orfLength}
\title{Length of ORF}
\usage{
orfLength(orf.table, aa = FALSE)
}
\arguments{
\item{orf.table}{A GFF-formatted \code{tibble}.}

\item{aa}{Logical, length in amino acids instead of bases.}
}
\value{
A vector of lengths.
}
\description{
Computing the lengths of all ORFs in an \code{orf.table}.
}
\details{
By default, computes the length of an ORF in bases, including the
stop codon. However, if \code{aa = TRUE}, then the length is in amino acids
after translation. This aa-length is the base-length divided by 3 and minus 
1, unless the ORF is truncated and lacks a stop codon.
}
\examples{
# See the example in the Help-file for findOrfs.

}
\seealso{
\code{\link{findOrfs}}.
}
\author{
Lars Snipen and Kristian Hovde Liland.
}
