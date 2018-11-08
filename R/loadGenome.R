#' loadGenome
#'
#' Loads BSgenome sequences for working with genomic coordinates
#'
#' @param genome.ver Reference genome (hg18, hg19, hg38)
#' @param exclude.scaffold.chrs If true (default), do not include scaffold chromosomes in the returned DNAStringSet object
#' @return DNAStringSet object with the chosen reference genome associated to genome.ver.
#' By default, the sequence of the '+' strand is downloaded.
#' @keywords Reference genome
#' @export
#' @examples
#' genome <- loadGenome (genome.ver = "hg38", exclude.scaffold.chrs = TRUE)
#'
loadGenome <- function (genome.ver="hg38", exclude.scaffold.chrs=TRUE){
  chrs <- paste0("chr", c(1:22, "X", "Y"))
  # Load the selected version of the genome
  if (genome.ver == "hg38")	{
    genome <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, chrs, strand = '+')
  } else if (genome.ver == "hg19"){
    genome <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chrs, strand=  '+')
  }
  else if (genome.ver == "hg18"){
    genome <- Biostrings::getSeq(BSgenome.Hsapiens.UCSC.hg18::Hsapiens, chrs, strand=  '+')
  }else stop("Specified genome.ver is not supported")

  # if (exclude.scaffold.chrs)
  #   genome <- genome[!grepl("_", names(genome))]

  return (genome)
}
