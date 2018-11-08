#' loadSignature
#'
#' Loads a signature of genomic coordinates and convert to destination genome reference version (if required)
#'
#' @param file Name of file containing signature coordinates (cvs format)
#' @param numberFeatures Limit the features to import from file to the first numberFeatures features;
#' if NULL (default) all the features in file are imported
#' @param coord.origin String with the reference genome version of the coordinates in file (supports values "hg38", "hg19", "hg18")
#' @param coord.dest If different than coord.origin, the coordinates in file will be converted from coord.origin to coord.dest.
#' Current support only for conversions from "hg38" to "hg19" and "hg18".
#' @return  GRanges object with the coordinates of the signature (extra columns are added as metadata to the GRanges object)
#' @keywords signatures
#' @export
#' @examples
#' signature <- loadSignature ()

loadSignature <- function (file, numberFeatures=0, coord.origin = "hg38", coord.dest="hg38" ){

  signature.orig <- GenomicRanges::makeGRangesFromDataFrame(read.csv(file))
  if (numberFeatures>0) {
  signature.orig <- signature.orig[1:numberFeatures]
}
  if (coord.origin == coord.dest)
    signature.dest <- signature.orig

  else{

    # need to convert biomarkers coordinates
    # Current support only for conversions from hg.38 to hg.19-hg.18
    ch.hg38.to.hg19 <- rtracklayer::import.chain("/media/mendel/reference_genomes/liftOver/hg38ToHg19.over.chain")
    signature.hg19 <- rtracklayer::liftOver(signature.orig, ch.hg38.to.hg19)	#returns a GRangesList object, since the mapping might be one-to-many
    #hist(elementNROWS(signature.hg19))								#check cardinality of the mapping (1-to-0 and 1-to-1 are the most popular)
    signature.hg19 <- unlist(signature.hg19)						#remove nesting, each coordinate is a new biomarker (potential overrepresentation of a hg38 biomarker that has one-to-many mappings to hg19)

    if (coord.dest =="hg19"){
      signature.dest <- signature.hg19
    }else{
      ch.hg19.to.hg18 <- rtracklayer::import.chain("/media/mendel/reference_genomes/liftOver/hg19ToHg18.over.chain")
      signature.hg18 <- rtracklayer::liftOver(signature.hg19, ch.hg19.to.hg18)
      signature.hg18 <- unlist(signature.hg18)						#remove nesting, each coordinate is a new biomarker (potential overrepresentation of a hg38 biomarker that has one-to-many mappings to hg19)
      signature.dest <- signature.hg18
    }

    return (signature.dest)
  }

}
