% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fastq.R
\name{readFastq}
\alias{readFastq}
\alias{writeFastq}
\title{Read and write FASTQ files}
\usage{
readFastq(in.file)
writeFastq(fdta, out.file)
}
\arguments{
\item{in.file}{url/directory/name of (gzipped) FASTQ file to read.}

\item{fdta}{FASTQ object to write.}

\item{out.file}{url/directory/name of (gzipped) FASTQ file to write.}
}
\value{
\code{\link{readFastq}} returns a \code{\link{tibble}} with the contents of the (gzipped) FASTQ
file stored in three columns of text. The first, named \samp{Header}, contains
the headerlines, the second, named \samp{Sequence}, contains the sequences and the third, named 
\samp{Quality} contains the base quality scores.

\code{\link{writeFastq}} produces a (gzipped) FASTQ file.
}
\description{
Reads and writes files in the FASTQ format.
}
\details{
These functions handle input/output of sequences in the commonly used FASTQ format,
typically used for storing DNA sequences (reads) after sequencing. If
filenames (\code{in.file} or \code{out.file}) have the extension \code{.gz} they will automatically be
compressed/uncompressed. 

The sequences are stored in a \code{\link{tibble}}, opening up all the possibilities in R for
fast and easy manipulations. The content of the file is stored as three columns, \samp{Header},
\samp{Sequence} and \samp{Quality}. If other columns are added, these will be ignored by
\code{\link{writeFastq}}.
}
\note{
These functions will only handle files where each entry spans one single line, i.e. not the
(uncommon) multiline FASTQ format.
}
\examples{
\dontrun{
# We need a FASTQ-file to read, here is one example file:
fq.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fastq.gz")

# Read and write
fdta <- readFastq(fq.file)
ok <- writeFastq(fdta[1:3,], out.file = "delete_me.fq")

# Make use of dplyr to copy parts of the file to another file
readFastq(fq.file) \%>\% 
  mutate(Length = str_length(Sequence)) \%>\% 
  filter(Length > 200) \%>\% 
  writeFasta(out.file = "long_reads.fa") # writing to FASTA file
}

}
\seealso{
code{\link{readFasta}}.
}
\author{
Lars Snipen and Kristian Hovde Liland.
}
\keyword{FASTQ}
\keyword{sequence}
