#' Searches for matches of a nucleotide sequence genome wide and returns matching coordinates
#'
#' @param genome DNAStringSet object with the reference genome (typically as provided by getSeq from package BSgenome)
#' @param pattern DNA sequence to search for matches.
#' @return GRanges object with the coordinates of every match for the pattern sequence in the genome (genome version as provided in genome.ver)
#' @keywords CG GRanges
#' @export
#' @examples sequenceCoordinates("hg38", "CG")  #searches for dinucleotides CpG genome-wide in hg38.

sequenceCoordinates <- function(genome, pattern ="CG"){

  matches <- Biostrings::vmatchPattern (Biostrings::DNAString (pattern), genome)
  matches <- GenomicRanges::GRanges (matches)
  return (matches)
}
