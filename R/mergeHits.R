#' Gnerate GRanges object with number of signatures, number of dinucleotides CG
#' number of N and N proportion per bin
#' @param gr.hits GRanges object with number of hits per bin as metadata
#' @param gr_CG.hits GRanges object with number of dinucleotides CG per bin as metadata
#' @param gr_N.hits Granges object with number of N per bin as metadata
#' @return GRanges object
#' @keywords GRanges
#' @export
#' @examples mergeHits(gr.hits, gr_CG.hits, gr_N.hits)

mergeHits <- function(gr.hits, gr_CG.hits, gr_N.hits){
  Nprop <- data.frame(GenomicRanges::elementMetadata(gr_N.hits)[["N"]]/GenomicRanges::width(GITAR.GSM.hg38.N.hits))
  gr_CG.hits$bin.type <- NULL
  gr_N.hits$bin.type <- NULL
  colnames(Nprop)<- "N.prop"
  GenomicRanges::mcols(gr.hits)=cbind(elementMetadata(gr.hits), elementMetadata(gr_CG.hits),
                                                elementMetadata(gr_N.hits), Nprop)
  return(gr.hits)



}
