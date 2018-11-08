#' Generate a Granges object (only seqname, start and end) to a BED-like file (3 fields format)
#'
#' @param gr GRanges object
#' @param file Name of BED file to be generated
#' @param chr.header Name of header field for chromosome
#' @param start.header Name of header field for start coordinate
#' @param end.header Name of Header field for end coordinate
#' @return BED file
#' @keywords BED GRanges
#' @export
#' @examples
#' Granges2Bed(gr, 'my_bed_file.bed')


GRanges2bed <- function (gr, file, chr.header = "chr", start.header = "start",
                         end.header = "end"){
  df <- data.frame(seqnames=seqnames(gr),
                   starts=start(gr)-1,
                   ends=end(gr))

  write.table(df, file=file, quote=F, sep="\t", row.names=F, col.names=F)

}
