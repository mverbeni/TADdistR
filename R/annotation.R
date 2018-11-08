#' Get Gencode GTF data from RAW GENCODE file
#' @param file File with genome annotation
#' @param genes.only If TRUE (default) only genes are considered
#' @return Annotation
#' @export

getGencodeGTF <- function(file, genes.only=TRUE) {
  contents <- rtracklayer::import(file)
  if (genes.only) {
    contents <- contents[contents$type == "gene"]
  }
  return(contents)
}

#' annotateBins
#'
#' This function loads two GRanges objects, one with a partition of the genome (disjoint GRanges object) and one annotation track with annotation
#'  of interest (e.g. genes, lncRNAs)
#' The function adds metadata to the bins GRanges object with the annotations falling in each bin of the partition.
#'
#' @param bins GRanges object with a disjoint partition of the genome
#' @param annotation GRanges object with annotations (e.g. genes, lncRNAs)
#' @param column.name Name of the field from annotation to be added as metadata column to the result GRanges.
#' @param ignore.strand Ignore strandedness for findOverlaps
#' @return GRanges object with a copy of bins with an additional metadata column named "overlaps" with the annotations from annotation falling into each bin
#' @keywords GRanges hits binning
#' @export
#' @examples
#' annotateBins(bins, annotation)

annotateBins <- function(bins, annotation, id.name="id", column.name="coding_genes", ignore.strand=TRUE){
  olaps <- GenomicRanges::findOverlaps(bins, annotation, ignore.strand=ignore.strand, minoverlap=1L, maxgap=-1L)
  f1 <- factor(S4Vectors::queryHits(olaps), levels=seq_len(S4Vectors::queryLength(olaps)))
  s <- splitAsList(mcols(annotation)[[id.name]][subjectHits(olaps)], f1)
  nhits <- sapply(s, length)
  s.asstring <- sapply(s, paste, collapse=",")
  mcols(bins)[[paste("n", column.name, sep=".")]] <- nhits
  mcols(bins)[[column.name]] <- s.asstring

  bins
}

#' Annotate bins per category
#' @param bins GRanges object with a disjoint partition of the genome
#' @param annotations GRanges object with annotations (e.g. genes, lncRNAs)
#' @return Annotation according to category
#' @export
annotateBinsFromGRanges <- function(bins, annotations) {

  bins.annotated <- annotateBins(bins, annotations[annotations$gene_type == "protein_coding"], id.name="gene_name", column.name="protein_coding")
  bins.annotated <- annotateBins(bins.annotated, annotations[annotations$gene_type == "lincRNA"], id.name="gene_name", column.name="lincRNA")
  bins.annotated <- annotateBins(bins.annotated, annotations[annotations$gene_type == "pseudogene"], id.name="gene_name", column.name="pseudogene")
  bins.annotated <- annotateBins(bins.annotated, annotations[annotations$gene_type == "antisense"], id.name="gene_name", column.name="antisense")
  bins.annotated <- annotateBins(bins.annotated, annotations[annotations$gene_type == "miRNA"], id.name="gene_name", column.name="miRNA")

  return(bins.annotated)
}

