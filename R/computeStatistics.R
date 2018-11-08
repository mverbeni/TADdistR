#' ComputeProbability
#'
#' Compute observed, expected and probability for a GRanges object given a counting for hits (features) and CGs.
#' Find bins with a number of signatuire significatly different from the expected,
#' calculated as N*p_0, where N is the total nuymber of bins, and the probability
#' p_0 is given by the number of CG in each bin divided by the total number
#' of CG in the reference genome
#'
#' @param gr GRanges object with a disjoint partition of the genome and
#'  two metadata columns named as specified in hist.colum.name and CG.column.name
#' @param hits.column.name Number of hits in each bin
#' @param CG.column.name Number of dinucleotides CG per bin
#' @param N.column.name Number of N per bin
#' @param N.prop.name Proportion of N in each bin
#' @param thrN Threshold for the fraction of N in each bin
#' @return Copy of GRanges object gr with new metadata columns "observed" "probability" and "expected" (previous metadata columns in gr are also kept)
#' @keywords GRanges probability observations expected
#' @export
#' @examples
#' computeStatistics (gr)

computeStatistics <- function (gr, hits.column.name="Lynch.wilcox.0.005.hits",
                               CG.column.name="CG.hits",N.column.name="N.hits",N.rate.column.name="N.prop", thrN = 0.0){

  numCG0 <- gr$CG.hits
  remove_cg0 <- which(numCG0==0)
  gr <- gr[-remove_cg0]
  if (thrN > 0) {
    N.rate <- as.numeric(elementMetadata(gr)[["N.rate"]])
    remove_n <- which(N.rate >=thrN)
    gr <- gr[-remove_n]
  }
  numCG <- gr$CG.hits
  nbins <- length(gr$Lynch.wilcox.0.005.hits)
  totCG <- sum(numCG)
  prob0 <- numCG/totCG
  hits <- gr$Lynch.wilcox.0.005.hits
  nhits <- sum(hits)
  nbins <- length(hits)
  expt <- nhits*prob0

  btest <- function(x,n,p){
    test <- binom.test(x,n,p, alternative = "greater")
    pv <- test$p.value
    return(pv)
  }
  resbtest <- rep(0,length(hits))
  #sl <- 0.05
  for (i in 1:nbins){
    pvalue <- btest(hits[i],nhits,prob0[i])
    # if (pvalue < sl) {
    resbtest[i] <- pvalue
    # }
  }
  norm.width <- width(gr)/100000
  hits.per.100K <- hits/norm.width
  dataCG_hits.df <- cbind.data.frame(hits,prob0,expt, norm.width, hits.per.100K,resbtest)
  colnames(dataCG_hits.df) <- c("observed", "probability","expected", "norm.width.100Kb", "hits.per.100Kb","p value")
  mcols(gr)=cbind(mcols(gr), dataCG_hits.df)
  grord <- gr[order(gr$'p value'),]
  return (grord)
}
