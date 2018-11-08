#' Generate a GRanges object from a BED file
#'
#' This function loads a BED-like file and stores it as a GRanges object.
#' The tab-delimited file must not be ordered, but must include a header line indicating the fields which includes chr.field, start.field,
#' end.field (minimum), and optionally: strand.field.
#' Any columns apart from the four columns described are loaded into the GRanges object as metadata.
#'
#' @param file Name of the BED file
#' @param chr.field Name in the header of the BED file associated with the chr field
#'  for the GRanges conversion
#' @param start.field Name in the header of the BED file associated with the start field
#' for the GRanges conversion
#' @param end.field Name in the header of the BED file associated with the end field
#' for the GRanges conversion
#' @param strand.field Name in the header of the BED file associated with the strand field
#' for the GRanges conversion
#' @return GRanges object
#' @keywords BED GRanges
#' @export
#' @examples
#' bed2Granges('my_bed_file.bed')

bed2GRanges <- function(file, chr.field = "chr", start.field = "start",
                        end.field = "end", strand.field = "strand"){
  df <- read.delim(file, header=T)

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE,
                                 seqnames.field = c(chr.field),
                                 start.field = start.field,
                                 end.field = end.field,
                                 strand.field = strand.field)
  return (gr)
}





