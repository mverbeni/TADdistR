#' TADbinning
#'
#' This function loads a BED-like file and stores it into a GRanges object. This BED-like file contains bin (e.g TAD) coordinates.
#' The function prints basic statistics on the width and number of bins the genome is divided into.
#' The tab-delimited file must not be ordered, but must include a header line indicating the fields which includes chr.field, start.field,
#' end.field (minimum), and optionally: strand.field.
#' Any columns apart from the four columns described are loaded into the GRanges object as metadata.
#'
#' @param file Name of BED file containing TAD coordinates
#' @param genome DNAStringSet object with the reference genome as returned by the function loadGenome
#' @param chr.field Field name in header associated to the chr field for the GRanges conversion
#' @param start.field Field name in header associated to the start field for the GRanges conversion
#' @param end.field Field name in header associated to the end field for the GRanges conversion
#' @param strand.field Field name in header associated to the strand field for the GRanges conversion
#' @param gap.consecutive.bins Size of the gap between two consecutive bins. I.e. if there are two bins a: [a1,a2], b: [b1,b2],
#' what is the value of b1-a2, so we consider the bins consecutive? (default 0). Set gap.consecutive.bins to 0 when consecutive gaps
#' share end/start coordinate (i.e. b1==a2))
#' @param verbose TRUE (default) plots and shows statistics on distributions of bins width
#' @return GRanges object with a whole genome partition/binning defined by the provided TADs coordinates
#' @keywords BED GRanges TAD binning
#' @export
#' @examples TADbinning('HiCtool_allchrs_TADS.bed', genome, chr.field="chr", start.field="start", end.field="end", gap.consecutive.bins=0, verbose=F)

TADbinning <- function(file, genome, chr.field="chr", start.field="start", end.field="end", strand.field="strand", gap.consecutive.bins=0, verbose=TRUE){

  # 1- Load to GRanges, sort, find gaps, merge overlapping domains and list domains and gaps as bins

  TADs.Granges <- bed2GRanges (file)

  if (gap.consecutive.bins == 0)			#consecutive bins share the end/start coordinates. Reduce all ranges' end limits in 1 to prevent these bins from being merged
    GenomicRanges::end(TADs.Granges) <- GenomicRanges::end(TADs.Granges)-1

  TADs.Granges.sorted <- GenomicRanges::sort(TADs.Granges)
  TADs.Granges.gaps <- GenomicRanges::gaps (TADs.Granges.sorted) # to eliminate this, need to change below at point 2-
  TADs.Granges.sorted.reduced <- GenomicRanges::reduce(TADs.Granges.sorted, min.gapwidth=0)		#argument min.gapwith set to 0 to prevent consecutive bins from being merged together
  TADs.Granges.reduced.gaps <- GenomicRanges::gaps (TADs.Granges.sorted.reduced)

  #add one last gap to the end of each chromosome
  TADs.Granges.ending.gaps <- GenomicRanges::GRanges()
  for (i in GenomeInfoDb::seqlevels(TADs.Granges.sorted.reduced)){
    #last bin:
    last <- tail(TADs.Granges.sorted.reduced[GenomeInfoDb::seqnames(TADs.Granges.sorted.reduced)==i], n=1)
    #new bin from the end of the last one to the end of the chromosome:
    TADs.Granges.ending.gaps <- suppressWarnings(append(TADs.Granges.ending.gaps,
                                      GenomicRanges::GRanges (seqnames=i, ranges=IRanges::IRanges(start=BSgenome::end(last)+1, end = BSgenome::width(genome[i]))) ))
  }

  TADs.Granges.sorted.reduced$bin.type <- rep("TAD", length(TADs.Granges.sorted.reduced))
  TADs.Granges.reduced.gaps$bin.type <- rep("GAP", length(TADs.Granges.reduced.gaps))
  TADs.Granges.ending.gaps$bin.type <- rep("GAP", length(TADs.Granges.ending.gaps))
  partition <- GenomicRanges::sort(c(TADs.Granges.sorted.reduced, TADs.Granges.reduced.gaps, TADs.Granges.ending.gaps))


  # 2- Plot and show statistics

  if (verbose){
    par(mfrow=c(3, 2))

    print ('Computing domain and gap sizes for ')
    print (file)

    TADs.Granges.sorted.width <- width(TADs.Granges.sorted)
    hist(TADs.Granges.sorted.width[TADs.Granges.sorted.width<(mean(TADs.Granges.sorted.width)*5)],
         main="Width of domains",
         xlab="Width of domain (bp)", ylab="Frequency")
    print ("Basic statistics on the width of genomic ranges")
    print(summary(TADs.Granges.sorted))
    print (summary (TADs.Granges.sorted.width))

    TADs.Granges.gaps.width <- width(TADs.Granges.gaps)
    hist(TADs.Granges.gaps.width[TADs.Granges.gaps.width<(mean(TADs.Granges.gaps.width)*5)],
         main="Distance between domains",
         xlab="Distance between domains (bp)", ylab="Frequency")
    print ("Basic statistics on the width of the SPANS in between every genomic range")
    print(summary(TADs.Granges.gaps))
    print (summary(TADs.Granges.gaps.width))

    TADs.Granges.sorted.reduced.width <- width(TADs.Granges.sorted.reduced)
    hist(TADs.Granges.sorted.reduced.width[TADs.Granges.sorted.reduced.width<(mean(TADs.Granges.sorted.reduced.width)*5)],
         main="Width of domains (combined)",
         xlab="Width of domain (bp)", ylab="Frequency")
    print ("Basic statistics on the width of COMBINED overlapping genomic ranges")
    print(summary(TADs.Granges.sorted.reduced))
    print (summary (TADs.Granges.sorted.reduced.width))

    TADs.Granges.reduced.gaps.width <- width(TADs.Granges.reduced.gaps)
    hist(TADs.Granges.reduced.gaps.width[TADs.Granges.reduced.gaps.width<(mean(TADs.Granges.reduced.gaps.width)*5)],
         main="Distance between domains (combined)",
         xlab="Distance between domains (bp)", ylab="Frequency")
    print ("Basic statistics on the width of the SPANS in between every COMBINED overlapping genomic range")
    print(summary(TADs.Granges.reduced.gaps))
    print(summary(TADs.Granges.reduced.gaps.width))

    partition.width <- width(partition)
    hist(partition.width[partition.width<(mean(partition.width)*5)],
         main="Bin size",
         xlab="Bin size (bp)", ylab="Frequency")
    print ("Bin size")
    print(summary(partition))
    print(summary(partition.width))
  }

  return (partition)

}
