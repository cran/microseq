#' @name readFastq
#' @title Read FASTQ files
#' @aliases readFastq
#' 
#' @description Reads biological sequences (DNA) in the FASTQ format.
#' 
#' @usage readFastq(in.file)
#' 
#' @param in.file url/directory/name of FASTQ file to read.
#' 
#' @details The FASTQ format is typically used for storing DNA sequences (reads) together
#' with quality scores from a sequencer. Reading a FASTQ-file with this function will result
#' in ignoring the quality score string, and the sequences are stored in a \code{Fasta}
#' object.
#' 
#' This function uses the \code{\link{readBStringSet}} function of the \code{Biostrings} package.
#' 
#' @return \code{\link{readFastq}} returns a \code{Fasta} object with the header and sequence part of the 
#' FASTQ file.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{readFasta}}, \code{\link{plot.Fasta}}, \code{\link{summary.Fasta}}.
#' 
#' @examples 
#' \dontrun{
#' ex.file <- file.path(file.path(path.package("microseq"),"extdata"),"small.fastq")
#' fdta <- readFastq(ex.file)
#' summary(fdta)
#' }
#' 
#' @keywords sequence FASTQ Fastq
#' 
#' @importFrom Biostrings readBStringSet
#' 
#' @export readFastq
#' 
readFastq <- function( in.file ){
  xss <- readBStringSet( filepath=in.file, format="fastq" )
  seq <- as.character( xss )
  names( seq ) <- NULL
  
  fdta <- data.frame(
    Header   = names( xss ),
    Sequence = seq,
    stringsAsFactors=F )
  class( fdta ) <- c( "Fasta", "data.frame" )
  return( fdta )
}



