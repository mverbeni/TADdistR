#' hitsPerBin
#'
#' This function loads two GRanges objects, one with a partition of the genome (disjoint GRanges object) and one with loci of interest (e.g. relevant biomarkers).
#' The function counts how many of the loci of interest (biomarkers) fall in each bin of the partition.
#'
#' @param bins GRanges object with a disjoint partition of the genome
#' @param features GRanges object with a list of loci of interest (e.g. genetic markers)
#' @param column.name Name of the metadata column to be added to the resulting GRanges with the number
#'  of features falling into each bin
#' @return GRanges object with a copy of bins including an additional metadata column 'hits'
#'  with the number features ranges falling into each bin.
#' @keywords GRanges hits binning
#' @export
#' @examples hitsPerBin(subject, query)

hitsPerBin <- function(bins, features, column.name = "hits"){
  hits<- suppressWarnings(GenomicRanges::findOverlaps (features, bins, minoverlap=1L, maxgap=-1L))
  GenomicRanges::mcols(bins)[[column.name]] <- S4Vectors::countRnodeHits(hits)
  return (bins)
}
